
set_plot, 'ps'
device, filename='sfr_bouwens.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 0, 0, 255, 4
TvLCT, 255, 128, 0, 5
!P.charsize = 1.2
!P.charthick = 1
!P.thick = 1 

; Plot Hopkins and Beacom data points. 
readcol, '../data/sfrsfr.data', rs, rserr, sf, sferr, /silent
sfe = 10.0^sf 
plotsym, 0, 0.5, /FILL
plot, rs, sfe, psym=8, /ylog, /xlog, xrange=[5.0e-1,30], $
      ytitle='!6SFR (M!D!9n!X!N yr!E-1!N Mpc!E-3 !N)', $
      xstyle=1, yrange=[1.0e-4,1], ytickformat='Exponent', $
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
oplot, rs_bouwens, sfr_bouwens, psym=8, color=4 
dy = sfr_bouwens-10.0^(y+dym)
dy2 = 10.0^(y+dyp)-sfr_bouwens 
dx2 = dxp 
dx1 = -dxm 
plot_err, rs_bouwens, sfr_bouwens, dy, dy2=dy2, dx1=dx1, dx2=dx2, color=4
x0 = rs_bouwens[0]
y0 = sfr_bouwens[0] 
x1 = rs_bouwens[0]
y1 = 1.5e-4
arrow, x0, y0, x1, y1, /data, hsize=200.0, color=4

; Plot our result. 
restore, 'sfrfiletemplate.sav'
sfrdata = read_ascii('set10/sfr.out', template=sfrfiletemplate)
sfr_tot = sfrdata.total_sfr     ; Msun yr^-1 Mpc^-3 
z = sfrdata.redshift
;oplot, z, sfr_tot 

; Add legend.
legend, ['Hopkins and Beacom 06','Bouwens et al. 11 (>0.06L!D*,z=3!N)'], $
        linestyle=[0,0], psym=[8,8], color=[2,4], charsize=1, /bottom

; Produce our data in alternative way.
limit = 160

restore, 'halo_template.sav'
n_data = read_ascii('set10/nofmh.out', template=stars_template)
rs_alt = z
sfr_alt = sfr_tot
nm = n_data.field001[0:*,0]
nh = size(nm, /n_elements)
for i = 0, 99 do begin 
   sfr_alt[i]=0.0
   nm = n_data.field001[0:*,i]
   for j = 1, nh-1 do begin 
      sfr_alt[i] = sfr_alt[i] + nm[j]
   endfor
   sfr_alt[i] = sfr_alt[i]*1.0e10 
endfor
oplot, rs_alt, sfr_alt

restore, 'halo_template.sav'
n_data = read_ascii('set10/nofmh.out', template=stars_template)
rs_alt = z
sfr_alt = sfr_tot
nm = n_data.field001[0:*,0]
nh = size(nm, /n_elements)
for i = 0, 99 do begin 
   sfr_alt[i]=0.0
   nm = n_data.field001[0:*,i]
   for j = limit, nh-1 do begin 
      sfr_alt[i] = sfr_alt[i] + nm[j]
   endfor
   sfr_alt[i] = sfr_alt[i]*1.0e10 
endfor
;oplot, rs_alt, sfr_alt, linestyle=2

readcol, 'set20/sfr.out', z, sfrtot, sfrp2, sfrp3, limsfr, format='f,d,d,d', /silent
limsfr = limsfr*3.0d3
oplot, z, limsfr, color=5, thick=4

device, /close_file
set_plot, 'X'

END 



