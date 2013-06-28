
; File: ldla.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Plots redshift evolution of l_DLA (i.e., dN/dX).  Reads data from
; one of the ldla*.dat files, which are created manually using the
; output of abr2.pro. 

set_plot, 'ps'
device, filename='ldla.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

;; window, xsize=1000, ysize=1000
;; Device, decomposed=0
!P.charsize = 1.5
!P.thick = 1.0

TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 0, 0, 255, 4

openr, lun, 'ldla4.dat', /get_lun
ldla_dat = fltarr(2,6)
readf, lun, ldla_dat
rs = ldla_dat[0,*]
ldla_avg = ldla_dat[1,*]
plotsym, 0, 1, /FILL
plot, rs, ldla_avg, xrange=[2,4.5], xtitle='!6z', ytitle='dN/dX', yrange=[0.0,0.15]
close, lun
free_lun, lun 

openr, lun, 'ldla_obs.dat', /get_lun 
ldla_obsdata = fltarr(6,6) 
readf, lun, ldla_obsdata
close, lun 
free_lun, lun 

rs = ldla_obsdata[2,*]
ldla_obs = ldla_obsdata[3,*]
oplot, rs, ldla_obs, psym=8, color=2

rserr = fltarr(2,6)
rserr[0,*] = rs[*] - ldla_obsdata[0,*] 
rserr[1,*] = ldla_obsdata[1,*] - rs[*]

ldla_err = fltarr(2,6) 
ldla_err[0,*] = ldla_obsdata[4,*]
ldla_err[1,*] = ldla_obsdata[5,*]

dy1 = ldla_err[1,*]
dy2 = ldla_err[0,*]

dx1 = rserr[0,*]
dx2 = rserr[1,*]

plot_err, rs, ldla_obs, dy1, dy2=dy2, dx1=dx1, dx2=dx2, color=2

readcol, '../data/ldla_noterdaeme.dat', zlow, zhigh, dz, dx, dndz, /silent 
z_noterdaeme = (zlow+zhigh)/2.0
dndx = dndz*(dz/dx)
oplot, z_noterdaeme, dndx, psym=8, color=4
dx1 = z_noterdaeme - zlow
dx2 = zhigh - z_noterdaeme 
dy = dx1*0.0 
dy2 = dx2 
dy2 = 0.0 
plot_err, z_noterdaeme, dndx, dy, dy2=dy2, dx1=dx1, dx2=dx2, color=4

legend, ['Noterdaeme et al. 2012', 'Prochaska et al. 2005'], linestyle=[0,0], color=[4,2], psym=[8,8]


device, /close_file
set_plot, 'X'

