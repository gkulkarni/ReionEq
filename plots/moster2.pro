
; File: moster2.pro
;  Cre: 2012
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $)

; Plots Moster relation.

;; set_plot, 'ps'
;; device, filename='moster.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0
window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4
!P.charsize = 2.0
!P.thick = 1.0
!P.charthick = 1

restore, 'halo_template.sav'
halo_data = read_ascii('set36/halos.out', template=stars_template)
stars_data = read_ascii('set36/halos_stars.out', template=stars_template)
halos = halo_data.field001[1:*,99]
stars = stars_data.field001[1:*,99]
mstar = stars/halos 

halos = halos*1.0e10 
stars = stars*1.0e10 
plot, halos, stars, /xlog, /ylog, yrange=[1.0e7,1.0e13],xrange=[1.0e10,1.0e15]
;xrange=[5.0e-1, 1.0e6], yrange=[1.0e-3,1.0e-1], xtitle='!NM!Dhalo!N/10!E10!N M!D' + sun + '!N', ytitle='!NM!D*!N/M!Dhalo!N'
;xyouts, 1.0e-1, 0.06, 'z = 0', charsize=2.0

gamma = 0.611
beta = 1.068
n = 0.02817
m1 = 10.0^11.9
moster = fltarr(270)
for i=0, 269 do begin moster[i]=2.0*n*halos[i]/((halos[i]/m1)^(-beta)+(halos[i]/m1)^gamma)
oplot, halos, moster, linestyle=2, color=2

; upper limit 
gamma = 0.611+0.012
beta = 1.068+0.051
n = 0.02817+0.00063
m1 = 10.0^(11.9+0.026)
moster = fltarr(270)
for i=0, 269 do begin moster[i]=2.0*n*halos[i]/((halos[i]/m1)^(-beta)+(halos[i]/m1)^gamma)
oplot, halos, moster, linestyle=2, color=2

; lower limit 
gamma = 0.611-0.01
beta = 1.068-0.044
n = 0.02817-0.00057
m1 = 10.0^(11.9-0.024)
moster = fltarr(270)
for i=0, 269 do begin moster[i]=2.0*n*halos[i]/((halos[i]/m1)^(-beta)+(halos[i]/m1)^gamma)
oplot, halos, moster, linestyle=2, color=2

