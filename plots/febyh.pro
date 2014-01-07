
PRO febyh, set, row 

  window, xsize=1000, ysize=1000
  Device, decomposed=0
  TvLCT, 255, 0, 0, 2 
  !P.charsize = 2

  fe_data2 = read_ascii(set + '/halos_fe.out', template=stars_template)
  gas_data2 = read_ascii(set + '/halos_gas.out', template=stars_template)
    
  lcolumn = 95 

  fe2 = fe_data2.field001[lcolumn:*,row]
  gas2 = gas_data2.field001[lcolumn:*,row]

  array_size = size(fe2)
  array_length = array_size[1] 

  febyh2 = fltarr(array_length)

  for i = 0, array_length-1 do begin 

     mH2 = gas2[i]*0.71
     
     if (fe2[i] eq 0.0) then begin 
        febyh2[i] = -10.0
     endif else begin 
        febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
     endelse

     print, i, febyh2[i] 

  endfor

  plot, febyh2 

END

