  window, xsize=1000, ysize=1000
  Device, decomposed=0
  TvLCT, 255, 0, 0, 2 
  !P.charsize = 2.0

  restore, 'halo_template.sav'
  tdyn_data = read_ascii('set60/halos_aux.out', template=stars_template)
  z = tdyn_data.field001[0,*]
  tdyn = tdyn_data.field001[250,*]
  plot, z, tdyn, /ylog, yrange=[1.0e5,1.0e10]

  tdyn = tdyn_data.field001[50,*]
  oplot, z, tdyn, thick=4

  readcol, 'tdyn.chk', z, tdyn2
  oplot, z, tdyn2

  r = tdyn2/tdyn 

