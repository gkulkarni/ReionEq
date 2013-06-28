
; File: reion.out 
;  Cre: 2012-09-14
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.3 $) 

; Can be used to plot reionization-related output of reion. Mainly
; gpi, but can be easily generalised. 

!P.multi = [0, 3, 0, 0, 1]
!P.charsize = 2.0
!P.thick = 1.0
Device, decomposed=0
TvLCT, 255, 0, 0, 2 

;; set_plot, 'ps'
;; device, filename='reion.ps', color=1, xsize=6.0, ysize=2.0, /Inches

; Plot gamma_pi 
restore, 'reionfiletemplate.sav'
reiondata = read_ascii('set32/reion.out', template=reionfiletemplate)
redshift = reiondata.z
avgtemp = reiondata.field07
; plot, redshift, avgtemp, /xlog, /ylog, xrange=[1,100], xtitle='!8z!6', ytitle='T (K)', charsize=1.2 
gpi = reiondata.field04
plot, redshift, gpi, /xlog, /ylog, xrange=[2,100], xstyle=1, xtitle='!6z', ytitle='log!D10!N(' + Greek('Gamma') + '!DPI!N/10!E-12!N s!E-1!N)' 

openr, lun, '../bh.dat', /get_lun
obsdata = fltarr(4,3)
readf, lun, obsdata
zobs = obsdata[0,*]
gobs = obsdata[1,*]
gerr_low = obsdata[2,*]
gerr_high = obsdata[3,*]
gerr_low[2] = 1.0e-14
gerr_low = gobs - gerr_low 
gerr_high[2] = gerr_high[2]+1.0e-13
gerr_high = gerr_high - gobs

oplot, zobs, gobs, psym=7, color=2
zerr = fltarr(3)
oploterror, zobs, gobs, gerr_high, errcolor=2, psym=3, /hibar
oploterror, zobs, gobs, gerr_low, errcolor=2, psym=3, /lobar

; Plot mass-metallicity relation at redshift two  
restore, 'halo_template.sav'
stars_data = read_ascii('set32/halos_stars.out', template=stars_template)
stars = stars_data.field001[1:*,95]
metal_data = read_ascii('set32/halos_obyh.out', template=stars_template)
metal = metal_data.field001[1:*,95]
stars = alog10(stars) - 1.0
metal = metal - 1.87 + 12.0 + 0.9
mets = smooth(metal, 30, /edge_truncate) 
sol = "156B
sun = '!9' + string(sol) + '!X'
plot, stars, mets, yrange=[7,12], xtitle='log!D10!N(M!D*!N/10!E10!N M!D' + sun + '!N)', ytitle='12+log!D10!N(O/H)'
xyouts, 0.0, 11.0, 'z=2'

;; stars_data = read_ascii('set32/halos_stars.out', template=stars_template)
;; stars = stars_data.field001[1:*,95]
;; metal_data = read_ascii('set32/halos_obyh.out', template=stars_template)
;; metal = metal_data.field001[1:*,95]
;; stars = alog10(stars) - 1.0
;; metal = metal - 1.87 + 12.0 + 0.4
;; oplot, stars, metal, color=3

openr, lun, '../data/erb_pointscorrect.dat', /get_lun
obsdata = fltarr(2,6)
readf, lun, obsdata
xobs = obsdata[0,*]
yobs = obsdata[1,*]
xobs = xobs-10.0
free_lun, lun
openr, lun, '../data/erb_xhigh.dat', /get_lun
xhigh = fltarr(2,6)
readf, lun, xhigh
free_lun, lun
openr, lun, '../data/erb_xlow.dat', /get_lun
xlow = fltarr(2,6)
readf, lun, xlow
free_lun, lun
a = xhigh[0,*]
b = xlow[0,*]
xerr = a-b
openr, lun, '../data/erb_yhigh.dat', /get_lun
yhigh = fltarr(2,6)
readf, lun, yhigh
free_lun, lun
openr, lun, '../data/erb_ylow.dat', /get_lun
ylow = fltarr(2,6)
readf, lun, ylow 
free_lun, lun
a = yhigh[1,*]
b = ylow[1,*]
yerr = a - b 
oploterror, xobs, yobs, xerr, yerr, errcolor=2, psym=3

; Plot Moster relation at redshift 0.
halo_data = read_ascii('set32/halos.out', template=stars_template)
halos = halo_data.field001[1:*,99]
stars = stars_data.field001[1:*,99]
mstar = stars/halos 
mss = smooth(mstar, 30, /edge_truncate) 
;plot, halos, mstar, /xlog, /ylog, xrange=[5.0e-1, 1.0e6],yrange=[1.0e-2,1.0e-1] 
plot, halos, mss, /xlog, /ylog, xrange=[5.0e-1, 1.0e6], yrange=[1.0e-2,1.0e-1], xtitle='!NM!Dhalo!N/10!E10!N M!D' + sun + '!N', ytitle='!NM!D*!N/M!Dhalo!N'
xyouts, 5.0e3, 0.06, 'z=0'

gamma = 0.6
beta = 1.5
n = 0.035
m1 = 10.0^11.5
moster = fltarr(270)
for i=0, 269 do begin moster[i]=2.0*n/((halos[i]*1.0e10/m1)^(-beta)+(halos[i]*1.0e10/m1)^gamma)
oplot, halos, moster, linestyle=2, color=2 

;; device, /Close_file 
;; set_plot, 'X'
