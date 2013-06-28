
PRO mzr_test4, halocolumn 

  window, xsize=1000, ysize=1000
  Device, decomposed=0
  TvLCT, 255, 0, 0, 2 
  TvLCT, 0, 127, 255, 3
  TvLCT, 255, 255, 0, 4
  !P.charsize = 2

  restore, 'halo_template.sav'
  stars_data = read_ascii('set49/halos_stars.out', template=stars_template)
  o_data = read_ascii('set49/halos_o.out', template=stars_template)
  gas_data =  read_ascii('set49/halos_gas.out', template=stars_template)

  hstrmass = fltarr(100)
  z = hstrmass
  o = hstrmass 

  for i = 0, 99 do begin 
     z[i] = stars_data.field001[0,i]
          
     hstrmass[i] = stars_data.field001[halocolumn,i]
     o[i] = o_data.field001[halocolumn,i]
  endfor 

  plot, z, hstrmass, /ylog, /xlog, xrange=[0.1,100], yrange=[1.0e-10,1.0e4], xstyle=1

  ;----------------------------------------------

  restore, 'halo_template.sav'
  stars_data = read_ascii('set10/halos_stars.out', template=stars_template)
  o_data = read_ascii('set10/halos_o.out', template=stars_template)
  gas_data =  read_ascii('set10/halos_gas.out', template=stars_template)

  hstrmass = fltarr(100)
  z = hstrmass

  for i = 0, 99 do begin 
     z[i] = stars_data.field001[0,i]

     hstrmass[i] = stars_data.field001[halocolumn,i]
     o[i] = o_data.field001[halocolumn,i]
  endfor 

  oplot, z, hstrmass, color=4

  ;----------------------------------------------

  restore, 'halo_template.sav'
  stars_data = read_ascii('set8/halos_stars.out', template=stars_template)
  o_data = read_ascii('set8/halos_o.out', template=stars_template)
  gas_data =  read_ascii('set8/halos_gas.out', template=stars_template)

  hstrmass = fltarr(100)
  z = hstrmass

  for i = 0, 99 do begin 
     z[i] = stars_data.field001[0,i]

     hstrmass[i] = stars_data.field001[halocolumn,i]
     o[i] = o_data.field001[halocolumn,i]
  endfor 

  oplot, z, hstrmass, color=3


END




