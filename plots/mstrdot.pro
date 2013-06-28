
PRO mstrdot, col 

  window, xsize=1000, ysize=1000
  Device, decomposed=0
  TvLCT, 255, 0, 0, 2 
  !P.charsize = 2

  restore, 'halo_template.sav'
  mstrdot_data = read_ascii('set41/halos_mstardot.out', template=stars_template)
  mstrdot = mstrdot_data.field001[col,*]*1.0e10 ; Msun / yr (Check!)
  z = mstrdot_data.field001[0,*] 
  plot, z, mstrdot, /xlog, /ylog, xrange=[3.0,100], yrange=[1.0e-6,1.0e1]

  vline, 7.0

  for i = 1, 20 do begin 
     mstrdot = mstrdot_data.field001[col+i*10,*]*1.0e10 ; Msun / yr (Check!)
     oplot, z, mstrdot
     if (i eq 5) then oplot, z, mstrdot, thick=5, color=2
  endfor

END 

 
