@ra_dec_cnvt

; program to extract the photometry for the Spitzer observations
; of the Webb calibrators

pro get_webbcal_phot,save_png=save_png,small_ap=small_ap

silent = 1
find_peak = 1

names = ['g191b2b','gd153','gd71','lds749b','hd209458','p041c','p177d','p330e', $
        'hd165459','1732526','1740346','1802271','1805292','1812095','1743045','snap2a','snap2b']
n_names = n_elements(names)

; bad regions to mask before photometry
bad_reg = intarr(n_names)
bad_reg_file = strarr(n_names)
bad_reg[0] = 1
bad_reg_file[0] = 'g191b2b_bad_ds9.reg'
;n_names = 1

; IRAC correction coefficients (PSF)
a = fltarr(6,4)
a[*,0] = [1.0114,-3.536E-6,-6.826E-5,-1.618E-8,1.215E-6,1.049E-6]
a[*,1] = [1.0138,8.401E-5,3.345E-7,1.885E-7,1.438E-6,1.337E-6]
a[*,2] = [1.0055,-3.870E-4,4.600E-5,1.956E-7,2.078E-6,9.970E-7]
a[*,3] = [1.0054,2.332E-4,-8.234E-5,-1.881E-7,6.520E-7,9.415E-7]

; IRAC correction coefficients (pixel phase)
irac_phase_a = [0.0535,0.0309]

if (not keyword_set(small_ap)) then begin
    irac_ap_rad = 10.
    irac_ap_sky = [12.,20.]
    ext_file_str = '_irac_rad10_sky12_20'
    irs_ap_rad = 10.
    irs_ap_sky = [12.,20.]
endif else begin
;    irac_ap_rad = 5.
;    irac_ap_sky = [10.,20.]
;    ext_file_str = '_irac_rad5_sky10_20'
    irac_ap_rad = 3.
    irac_ap_sky = [10.,20.]
    ext_file_str = '_irac_rad3_sky10_20'
    irs_ap_rad = 3.
    irs_ap_sky = [10.,20.]
endelse

; read in IRAC AOR ids
readcol,'irac_aorids.txt',xnames,xid1,xid2,xid3,xid4,format='(A,A,A,A,A)'
n_xnames = n_elements(xnames)
xid = strarr(4,n_xnames)
xid[0,*] = xid1
xid[1,*] = xid2
xid[2,*] = xid3
xid[3,*] = xid4

;irac_aorids = ['r28510976','r28511232','r28812544','r28812800','r28957696', $
;               'r7658240','r7655680','r7654912','r7656960','r11948800','r25130240']
;irac_indx = [1,5,2,3,0,4,0,11,10,6,7,8,9,0] - 1
irs_aorids = ['r28647424','r28647680','r28647936','r28648192','r28648448','r28648704','r14912256']
irs_indx = [0,0,0,0,7,6,0,0,0,1,2,3,4,5,0,0] - 1

;irac_indx = replicate(-1,n_names)
;irs_indx = replicate(-1,n_names)

; get coordinates

readcol,'coordinates.txt',cor_name,ra_str,dec_str,format='(A12,A12,A12)'
cor_name = strlowcase(cor_name)

openw,unit1,'webbcal_spitzer_phot'+ext_file_str+'.dat',/get_lun
printf,unit1,'#Spitzer photometry of Webb calibration stars'
printf,unit1,'#Karl Gordon, get_webbcal_phot.pro, 12 June 2009'
printf,unit1,'#'
printf,unit1,'#name, filter, xpos, ypos, flux[mJy], unc[mJy], total sky[mJy], total sky unc[mJy]'

openw,unit2,'webbcal_spitzer_indiv_phot'+ext_file_str+'.dat',/get_lun
printf,unit2,'#Spitzer photometry of Webb calibration stars'
printf,unit2,'#Derived from individual BCD images'
printf,unit2,'#Karl Gordon, get_webbcal_phot.pro, 12 June 2009'
printf,unit2,'#'
printf,unit2,'#name, filter, xpos, ypos, flux[mJy], unc[mJy], total sky[mJy], total sky unc[mJy] '

; get Spitzer data
;n_names = 1
;for k = 0,(n_names-1) do begin
for k = (n_names-2),(n_names-1) do begin
;for k = 5,5 do begin
    print,'working on ' + names[k]
    indxs = where(names[k] EQ cor_name,n_indxs)
    if (n_indxs EQ 0) then begin
        print,'no cor match for ' + names[k]
        stop
    endif
    ra_ref = 15.*ra_str2hr(ra_str[indxs[0]])
    dec_ref = dec_str2deg(dec_str[indxs[0]])

    ; see if there is IRAC data
    files = file_search('IRAC/data/*'+names[k]+'*.fits',count=n_files)
    if (n_files GT 0) then begin
        print,'found IRAC data n = '+strtrim(n_files,2)
        for i = 0,(n_files-1) do begin
;        for i = 3,(n_files-1) do begin
            fits_read,files[i],image,header
            print,files[i]
            save_image = image
            if (bad_reg[k]) then begin
                get_exclude_regions,'IRAC/data/'+bad_reg_file[k],exclude_regions,n_regions,silent=silent
                nan_dce_body,image,header,exclude_regions=exclude_regions,silent=silent
            endif
            extast,header,ast_info
            getrot,ast_info,rot,cdelt
            image_scale = 0.5*total(abs(cdelt)*3600.)
            ad2xy,ra_ref,dec_ref,ast_info,x_ref,y_ref
            if (keyword_set(save_png)) then make_png = repstr(files[i],'.fits','_phot') else make_png = 0
;            if (keyword_set(small_ap)) then begin ; find source
;                kgphot,save_image,x_ref,y_ref,obj_flux,obj_sn,sky_flux,obj_flux_un, $
;                       ap_rad=10.*image_scale,ap_sky=[10.,20.]*image_scale,make_png=make_png, $
;                       silent=silent,find_peak=find_peak,image_scale=image_scale
;            endif
            kgphot,image,x_ref,y_ref,obj_flux,obj_sn,sky_flux,obj_flux_unc,sky_flux_unc, $
                   ap_rad=irac_ap_rad*image_scale,ap_sky=irac_ap_sky*image_scale,make_png=make_png, $
                   silent=silent,find_peak=find_peak,image_scale=image_scale
            ch_pos = strpos(files[i],'ch')
            filter = 'IRAC'+strmid(files[i],ch_pos+2,1)
            printf,unit1,names[k],filter,x_ref,y_ref,obj_flux,obj_flux/obj_sn,sky_flux,sky_flux_unc, $
                   format='(A15,A6,F8.2,F8.2,E10.2,E10.2,E10.2,E10.2)'

            ; now get the photometry from the individual images
            indxs = where(xnames EQ names[k],n_indxs)
            if (n_indxs GE 0) then begin 
                chn_num = fix(strmid(files[i],ch_pos+2,1)) - 1
                ifiles = file_search('IRAC/'+xid[chn_num,indxs[0]]+'/ch'+strmid(files[i],ch_pos+2,1) + $
                                     '/bcd/*_bcd.fits',count=n_ifiles)
            endif else n_ifiles = 0

            if (n_ifiles GT 0) then begin
                ; get new ra,dec
                xy2ad,x_ref,y_ref,ast_info,new_ra,new_dec
                for m = 0,(n_ifiles-1) do begin
                    fits_read,ifiles[m],iimage,iheader
                    if (bad_reg[k]) then begin
                        get_exclude_regions,'IRAC/data/'+bad_reg_file[k],exclude_regions,n_regions,silent=silent
                        nan_dce_body,iimage,iheader,exclude_regions=exclude_regions,silent=silent
                    endif
                    iimage_size = size(iimage)
                    extast,iheader,iast_info
                    getrot,iast_info,irot,icdelt
                    iimage_scale = 0.5*total(abs(icdelt)*3600.)
                    ad2xy,new_ra,new_dec,iast_info,ix_ref,iy_ref
                    if ((ix_ref GE 0) AND (ix_ref LT iimage_size[1]) AND $
                        (iy_ref GE 0) AND (iy_ref LT iimage_size[2])) then begin
                        if (keyword_set(save_png2)) then make_png = repstr(files[i],'.fits','_phot'+strtrim(string(m+1),2)) $
                        else make_png = 0
                        kgphot,iimage,ix_ref,iy_ref,obj_flux,obj_sn,sky_flux,obj_flux_un, $
                               ap_rad=irac_ap_rad*iimage_scale,ap_sky=irac_ap_sky*iimage_scale,make_png=make_png, $
                               silent=silent,find_peak=0,image_scale=iimage_scale
                        ch_pos = strpos(files[i],'ch')
                        filter = 'IRAC'+strmid(files[i],ch_pos+2,1)
                        ; correction factor for distortions
                        corfac = a[0,chn_num] + a[1,chn_num]*(ix_ref - 128.) + a[2,chn_num]*(iy_ref - 128.) + $
                                 a[3,chn_num]*(ix_ref - 128.)*(iy_ref - 128.) + a[4,chn_num]*(ix_ref - 128.)^2 + $
                                 a[5,chn_num]*(iy_ref - 128.)^2
                        ; correction for pixel phase
                        if ((chn_num EQ 0) OR (chn_num EQ 1)) then begin
                            p_rad = sqrt(((ix_ref - fix(ix_ref)) - 0.5)^2 + ((iy_ref - fix(iy_ref)) - 0.5)^2)
                            corfac2 = 1.0 + irac_phase_a[chn_num]*((1./sqrt(2.*!PI)) - p_rad)
                        endif else corfac2 = 1.0
                        ; apply the correction
                        obj_flux *= corfac*corfac2
                        printf,unit2,names[k],m+1,filter,ix_ref,iy_ref,obj_flux,obj_flux/obj_sn,sky_flux,sky_flux_unc, $
                              format='(A15,I4,A6,F8.2,F8.2,E10.2,E10.2,E10.2,E10.2)'
                    endif
                endfor
            endif
        endfor
        
    endif

    ; see if there is IRS data
    files = file_search('IRS/data/*'+names[k]+'*.fits',count=n_files)
    if (n_files GT 0) then begin
        print,'found IRS data n = ' + strtrim(n_files,2)
        fits_read,files[0],image,header
        extast,header,ast_info
        getrot,ast_info,rot,cdelt
        image_scale = 0.5*total(abs(cdelt)*3600.)
        ad2xy,ra_ref,dec_ref,ast_info,x_ref,y_ref
        if (keyword_set(save_png)) then make_png = repstr(files[0],'.fits','_phot') else make_png = 0
        if (keyword_set(small_ap)) then begin ; find source
            kgphot,image,x_ref,y_ref,obj_flux,obj_sn,sky_flux,obj_flux_unc, $
                   ap_rad=10.*image_scale,ap_sky=[10.,20.]*image_scale,make_png=make_png, $
                   silent=silent,find_peak=find_peak,image_scale=image_scale
        endif
        kgphot,image,x_ref,y_ref,obj_flux,obj_sn,sky_flux,obj_flux_unc,sky_flux_unc, $
               ap_rad=irs_ap_rad*image_scale,ap_sky=irs_ap_sky*image_scale,make_png=make_png, $
               silent=silent,find_peak=find_peak,image_scale=image_scale
        filter = 'IRSB'
        printf,unit1,names[k],filter,x_ref,y_ref,obj_flux,obj_flux/obj_sn,sky_flux,sky_flux_unc, $
               format='(A15,A6,F8.2,F8.2,E10.2,E10.2,E10.2,E10.2)'

        ; now get the photometry from the individual images
        if (irs_indx[k] GE 0) then begin
            ifiles = file_search('IRS/'+irs_aorids[irs_indx[k]]+'/ch0/bcd/*_bcdb.fits', $
                                 count=n_ifiles)
        endif else n_ifiles = 0

            if (n_ifiles GT 0) then begin
                ; get new ra,dec
                xy2ad,x_ref,y_ref,ast_info,new_ra,new_dec
                for m = 0,(n_ifiles-1) do begin
                    fits_read,ifiles[m],iimage,iheader
                    iimage_size = size(iimage)
                    extast,iheader,iast_info
                    getrot,iast_info,irot,icdelt
                    iimage_scale = 0.5*total(abs(icdelt)*3600.)
                    ad2xy,new_ra,new_dec,iast_info,ix_ref,iy_ref
                    if ((ix_ref GE 0) AND (ix_ref LT iimage_size[1]) AND $
                        (iy_ref GE 0) AND (iy_ref LT iimage_size[2])) then begin
                        if (keyword_set(save_png2)) then make_png = repstr(files[0],'.fits','_phot'+strtrim(string(m+1),2)) $
                        else make_png = 0
                        kgphot,iimage,ix_ref,iy_ref,obj_flux,obj_sn,sky_flux,obj_flux_un,sky_flux_unc, $
                               ap_rad=irs_ap_rad*iimage_scale,ap_sky=irs_ap_sky*iimage_scale,make_png=make_png, $
                               silent=silent,find_peak=0,image_scale=iimage_scale
                        filter = 'IRSB'
                        printf,unit2,names[k],m+1,filter,ix_ref,iy_ref,obj_flux,obj_flux/obj_sn,sky_flux,sky_flux_unc, $
                              format='(A15,I4,A6,F8.2,F8.2,E10.2,E10.2,E10.2,E10.2)'
                    endif
                endfor
            endif
    endif
endfor
free_lun,unit1
free_lun,unit2

end
