
window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4 
TvLCT, 1, 3, 64, 5 
!P.charsize = 2
!P.multi=1

; Plot Hopkins and Beacom data points. 
readcol, '../data/sfrsfr.data', rs, rserr, sf, sferr, /silent
sfe = 10.0^sf 
plotsym, 0, 1, /FILL
plot, rs, sfe, psym=8, /ylog, /xlog, xrange=[2.0e-1,30], $
      ytitle='!6SFR (Msun yr!E-1!N Mpc!E-3 !N)', $
      xstyle=1, yrange=[1.0e-4,1], ytickformat='Exponent', $
      xtitle='redshift', /nodata 
oplot, rs, sfe, psym=8, color=2
yup = 10.0^(sf+sferr*0.5)
ylo = 10.0^(sf-sferr*0.5)
dy = yup-ylo
plot_err, rs, sfe, dy, dx1=rserr, color=2

; Plot our result. 
restore, 'sfrfiletemplate.sav'
sfrdata = read_ascii('set76/sfr.out', template=sfrfiletemplate)
sfr_tot = sfrdata.total_sfr ; Msun yr^-1 Mpc^-3 
z = sfrdata.redshift
oplot, z, sfr_tot 

sfrdata = read_ascii('set78/sfr.out', template=sfrfiletemplate)
sfr_tot = sfrdata.total_sfr ; Msun yr^-1 Mpc^-3 
z = sfrdata.redshift
oplot, z, sfr_tot, color=4

; Plot Tirth result. 
readcol, 'WMAP5atomic_sfr.dat', z, sfp2, sfp3, sft, rhot, /silent  
oplot, z, sfp2, linestyle=2

; Plot Hernquist and Springel 2003 result. 
sfr_hs = sfp2 
n = size(sfp2, /n_elements)
a = 0.012
b = 0.041
sfr_z0 = 0.013 ; Msun yr^-1 Mpc^-3 
omega_nr = 0.2646
omega_lambda = 0.734
for i = 0, n-1 do begin 
   x = (omega_nr*(1.0+z[i])^3+omega_lambda)^(1.0/3.0)
   sfr_hs[i] = sfr_z0 * x^2 / (1.0 + a*(x-1)^3*exp(b*x^(7.0/4.0))) ; Msun yr^-1 Mpc^-3 
endfor 
oplot, z, sfr_hs, linestyle=3

; Plot Hopkins and Beacom 2006 fit. 
sfr_hb = sfr_hs 
a = 0.0170
b = 0.13
c = 3.3
d = 5.3
smallh = 0.7 
for i = 0, n-1 do begin 
   sfr_hb[i] = (a+b*z[i])*smallh/(1.0+(z[i]/c)^d) ; Msun yr^-1 Mpc^-3 
endfor 
oplot, z, sfr_hb, linestyle=4 

; Add legend.
legend, ['my result', 'Tirth result', 'Hopkins and Beacom 06',$
         'Hernquist and Springel 03','Hopkins and Beacom 06 (fit)'], $
        linestyle=[0,2,0,3,4], psym=[0,0,8,0,0], color=[-1,-1,2,-1,-1], /bottom

END 



