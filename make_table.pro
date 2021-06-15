
pro make_table

file = 'webbcal_spitzer_ave_phot_all_rad3_sky10_20.dat'
readcol,file,format='(A15,A6,F8.2,F8.2,F10.2,F10.2,F10.2,F10.2,F8.2,I6)', $
  name,band,xpos,ypos,flux,flux_unc,sky,sky_unc,sn,num_images

openw,unit1,'webbcal_spitzer_phot.tex',/get_lun

sindxs = sort(name)
uindxs = uniq(name[sindxs])

n_uniq = n_elements(uindxs)
for i = 0,(n_uniq-1) do begin
    k = sindxs[uindxs[i]]
    indxs = where(name[k] EQ name,n_indxs)
    printf,unit1,format='(A12,$)',name[k] + ' &'

    bindx = where(band[indxs] EQ 'IRAC1',n_bindxs)
    if (n_bindxs EQ 1) then begin
        printf,unit1,format='(A,$)',string(flux[indxs[bindx[0]]],format='(E10.2)') + ' &'
        printf,unit1,format='(A,$)',string(flux_unc[indxs[bindx[0]]],format='(E10.2)') + ' &'
;        printf,unit1,format='(A,$)',string(num_images[indxs[bindx[0]]],format='(I4)') + ' &'
    endif else printf,unit1,format='(A,$)',' \nodata & \nodata &'
    
    bindx = where(band[indxs] EQ 'IRAC2',n_bindxs)
    if (n_bindxs EQ 1) then begin
        printf,unit1,format='(A,$)',string(flux[indxs[bindx[0]]],format='(E10.2)') + ' &'
        printf,unit1,format='(A,$)',string(flux_unc[indxs[bindx[0]]],format='(E10.2)') + ' &'
;        printf,unit1,format='(A,$)',string(num_images[indxs[bindx[0]]],format='(I4)') + ' &'
    endif else printf,unit1,format='(A,$)',' \nodata & \nodata &'
    
    bindx = where(band[indxs] EQ 'IRAC3',n_bindxs)
    if (n_bindxs EQ 1) then begin
        printf,unit1,format='(A,$)',string(flux[indxs[bindx[0]]],format='(E10.2)') + ' &'
        printf,unit1,format='(A,$)',string(flux_unc[indxs[bindx[0]]],format='(E10.2)') + ' &'
;        printf,unit1,format='(A,$)',string(num_images[indxs[bindx[0]]],format='(I4)') + ' &'
    endif else printf,unit1,format='(A,$)',' \nodata & \nodata &'
    
    bindx = where(band[indxs] EQ 'IRAC4',n_bindxs)
    if (n_bindxs EQ 1) then begin
        printf,unit1,format='(A,$)',string(flux[indxs[bindx[0]]],format='(E10.2)') + ' &'
        printf,unit1,format='(A,$)',string(flux_unc[indxs[bindx[0]]],format='(E10.2)') + ' &'
;        printf,unit1,format='(A,$)',string(num_images[indxs[bindx[0]]],format='(I4)') + ' &'
    endif else printf,unit1,format='(A,$)',' \nodata & \nodata &'
    
    bindx = where(band[indxs] EQ 'IRSB',n_bindxs)
    if (n_bindxs EQ 1) then begin
        printf,unit1,format='(A,$)',string(flux[indxs[bindx[0]]],format='(E10.2)') + ' &'
        printf,unit1,format='(A,$)',string(flux_unc[indxs[bindx[0]]],format='(E10.2)') + ' &'
;        printf,unit1,format='(A,$)',string(num_images[indxs[bindx[0]]],format='(I4)') + ' &'
    endif else printf,unit1,format='(A,$)',' \nodata & \nodata &'
;    endif else printf,unit1,format='(A,$)',' \nodata & \nodata & \nodata &'
    
    bindx = where(band[indxs] EQ 'MIPS1',n_bindxs)
    if (n_bindxs EQ 1) then begin
        printf,unit1,format='(A,$)',string(flux[indxs[bindx[0]]],format='(E10.2)') + ' &' 
        printf,unit1,format='(A)',string(flux_unc[indxs[bindx[0]]],format='(E10.2)') + ' \\'
;        printf,unit1,format='(A)',string(num_images[indxs[bindx[0]]],format='(I4)') + ' \\'
    endif else printf,unit1,format='(A)',' \nodata & \nodata \\'
    
    
endfor

free_lun,unit1

end
