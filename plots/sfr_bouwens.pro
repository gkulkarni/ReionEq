window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 255, 255, 0, 3
TvLCT, 255, 0, 255, 4
TvLCT, 0, 255, 255, 5
TvLCT, 0, 255, 0, 6
TvLCT, 201, 113, 0, 7

!P.charsize = 2.0
!P.charthick = 1
!P.thick = 1 

; Plot Hopkins and Beacom data points. 
readcol, '../data/sfrsfr.data', rs, rserr, sf, sferr, /silent
sfe = 10.0^sf 
plotsym, 0, 1, /FILL
plot, rs, sfe, psym=8, /ylog, /xlog, xrange=[5.0e-1,30], $
      ytitle='!6SFR (Msun yr!E-1!N Mpc!E-3 !N)', $
      xstyle=1, yrange=[1.0e-5,1], ytickformat='Exponent', $
      xtitle='z', /nodata 
oplot, rs, sfe, psym=8, color=2
yup = 10.0^(sf+sferr*0.5)
ylo = 10.0^(sf-sferr*0.5)
dy = yup-ylo
plot_err, rs, sfe, dy, dx1=rserr, color=2

; Plot Bouwens data points.
readcol, 'bouwens-data.dat', x, y, dxm, dxp, dym, dyp, /silent
rs_bouwens = x 
sfr_bouwens = 10.0^y ; Msun yr^-1 Mpc^-3 
oplot, rs_bouwens, sfr_bouwens, psym=8, color=5 
dy = sfr_bouwens-10.0^(y+dym)
dy2 = 10.0^(y+dyp)-sfr_bouwens 
dx2 = dxp 
dx1 = -dxm 
plot_err, rs_bouwens, sfr_bouwens, dy, dy2=dy2, dx1=dx1, dx2=dx2, color=5
x0 = rs_bouwens[0]
y0 = sfr_bouwens[0] 
x1 = rs_bouwens[0]
y1 = 1.5e-4
arrow, x0, y0, x1, y1, /data, hsize=13.0, color=5

; Plot our result. 
restore, 'sfrfiletemplate.sav'
sfrdata = read_ascii('set192/sfr.out', template=sfrfiletemplate)
sfr_tot = sfrdata.total_sfr*0.5     ; Msun yr^-1 Mpc^-3 
z = sfrdata.redshift
oplot, z, sfr_tot 

; Add legend.
legend, ['Hopkins and Beacom 06','Bouwens et al. 11 (>0.06L!D*,z=3!N)'], $
        linestyle=[0,0], psym=[8,8], color=[2,5], /bottom

;; Plot luminosity-limited SFR result.
readcol, 'set192/sfr.out', z, sfrtot, sfrp2, sfrp3, limsfr, format='f,d,d,d', /silent
oplot, z, limsfr, psym=-6, color=4

;; Plot luminosity-limited SFR result.
readcol, 'set192/sfr.out', z, sfrtot, sfrp2, sfrp3, limsfr, format='f,d,d,d', /silent
;oplot, z, limsfr, psym=-6, color=6

END 

