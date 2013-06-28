
PRO sfrt, set_name

  window, xsize=700, ysize=700
  Device, decomposed=0
  TvLCT, 255, 0, 0, 2 
  TvLCT, 0, 255, 0, 3
  !P.charsize=2 

  readcol, strtrim(set_name)+'/sfr.out', z, sfr, sfrp2, sfrp3, sfrlim 
  plot, z, sfr, /xlog, /ylog, xrange=[1,100] 
  oplot, z, sfrp2, thick=4, color=2
  oplot, z, sfrp3, thick=4, color=3

END

