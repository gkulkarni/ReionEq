
set_plot, 'ps'
device, filename='sfr_talk.ps', xsize=7.0, ysize=7.0, $
        /inches, color=1, /HELVETICA, yoffset=1.0
!P.font = 0 

TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
!P.charsize = 1.5
!P.thick = 3
!P.charthick = 1

restore, 'sfrfiletemplate.sav'
sfrdata = read_ascii('set68/sfr.out', template=sfrfiletemplate)
sfr_tot = sfrdata.total_sfr
z = sfrdata.redshift

readcol, '../data/sfrsfr.data', rs, rserr, sf, sferr 
sfe = 10.0^sf 

plotsym, 0, 1, /FILL
plot, rs, sfe, psym=8, /ylog, /xlog, xrange=[1,100], $
      ytitle='SFR (Msun yr!E-1!N Mpc!E-3 !N)', $
      xstyle=1, yrange=[1.0e-10,1], ytickformat='Exponent', $
      xtitle='redshift', xthick=3, ythick=3, /nodata 
oplot, rs, sfe, psym=8, color=2

yup = 10.0^(sf+sferr*0.5)
ylo = 10.0^(sf-sferr*0.5)
dy = yup-ylo
plot_err, rs, sfe, dy, dx1=rserr, color=2

oplot, z, sfr_tot 

legend, ['Hopkins and Beacom 2006'], psym=[8], color=[2], /bottom, box=0

device, /close_file
set_plot, 'X'

