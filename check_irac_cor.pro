
pro check_irac_cor,ix_ref,iy_ref

a = fltarr(6,4)
a[*,0] = [1.0114,-3.536E-6,-6.826E-5,-1.618E-8,1.215E-6,1.049E-6]
a[*,1] = [1.0138,8.401E-5,3.345E-7,1.885E-7,1.438E-6,1.337E-6]
a[*,2] = [1.0055,-3.870E-4,4.600E-5,1.956E-7,2.078E-6,9.970E-7]
a[*,3] = [1.0054,2.332E-4,-8.234E-5,-1.881E-7,6.520E-7,9.415E-7]

for chn_num = 0, 3 do begin
    corfac = a[0,chn_num] + a[1,chn_num]*(ix_ref - 128.) + a[2,chn_num]*(iy_ref - 128.) + $
      a[3,chn_num]*(ix_ref - 128.)*(iy_ref - 128.) + a[4,chn_num]*(ix_ref - 128.)^2 + $
      a[5,chn_num]*(iy_ref - 128.)^2
    print,chn_num+1, corfac
endfor

end

