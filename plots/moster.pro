
; File: moster.pro
;  Cre: 2012
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.2 $)

; Plots the Moster relation in our code.

PRO moster, set_name

  window, xsize=1000, ysize=1000
  Device, decomposed=0
  TvLCT, 255, 0, 0, 2 
  TvLCT, 0, 127, 255, 3
  TvLCT, 255, 255, 0, 4
  !P.charsize = 2.0

  restore, 'halo_template.sav'
  halo_data = read_ascii(strtrim(set_name)+'/halos.out', template=stars_template)
  stars_data = read_ascii(strtrim(set_name)+'/halos_stars.out', template=stars_template)
  halos = halo_data.field001[1:*,99]
  stars = stars_data.field001[1:*,99]
  mstar = stars/halos 
  mss = smooth(mstar, 30, /edge_truncate) 
  plot, halos, mss, /xlog, /ylog, xrange=[5.0e-1, 1.0e6], yrange=[1.0e-3,1.0e-1], $
        xtitle='!NM!Dhalo!N/10!E10!N M!D!9n!X!N', ytitle='!NM!D*!N/M!Dhalo!N'

  gamma = 0.6
  beta = 1.5
  n = 0.035
  m1 = 10.0^11.5
  moster = fltarr(270)
  for i=0, 269 do begin 
     moster[i]=2.0*n/((halos[i]*1.0e10/m1)^(-beta)+(halos[i]*1.0e10/m1)^gamma)
  endfor
  oplot, halos, moster, linestyle=2, color=2

END
