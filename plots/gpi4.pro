
PRO gpi4, set 

  set_plot, 'ps'
  device, filename='gpi_check2.ps', xsize=10.0, ysize=5.0, /inches, color=1, yoffset=1.0

  ;; window, xsize=2000, ysize=1000
  ;; Device, decomposed=0
  TvLCT, 255, 0, 0, 2 
  TvLCT, 0, 127, 255, 3
  TvLCT, 255, 255, 0, 4 
  TvLCT, 1, 3, 64, 5 
  TvLCT, 1, 3, 64, 6

  !P.charsize = 1
  !P.multi = [0,2,1,0,1]

  setnr = strtrim(uint(set),2)

; Plot SFR. 

; Plot Hopkins and Beacom data points. 
  readcol, '../data/sfrsfr.data', rs, rserr, sf, sferr, /silent
  sfe = 10.0^sf 
  plotsym, 0, 0.7, /FILL
  plot, rs, sfe, psym=8, /ylog, /xlog, xrange=[5.0e-1,30], $
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
  sfrdata = read_ascii('set'+setnr+'/sfr.out', template=sfrfiletemplate)
  sfr_tot = sfrdata.total_sfr   ; Msun yr^-1 Mpc^-3 
  z = sfrdata.redshift
  oplot, z, sfr_tot 

; Add legend.
  legend, ['Hopkins and Beacom 06'], linestyle=[0], psym=[8], color=[2], /bottom, charsize=1

;-------------------------------------------------------------------

; Plot Gamma_PI 

  restore, 'reionfiletemplate.sav'
  reiondata = read_ascii('set'+setnr+'/reion.out', template=reionfiletemplate)
  redshift = reiondata.z
  gpi = reiondata.field04
  gpi = gpi * 1.0e12 ;* 0.15
  ticks = [2,3,4,5,6,7,8,9,10,20]
  nticks = n_elements(ticks)
  plot, redshift, gpi, /ylog, xrange=[2,20], xstyle=1, xtitle='!6z', $
        ytitle='log!D10!N(!7C!6!DHI!N/10!E-12!Ns!E-1!N)', $
        yrange=[1.0e-4, 1.0e2], ytickformat='exponent', /xlog , $
        xticks=nticks-1, xtickv=ticks

  plotsym, 0, 0.7, /FILL
  zerr = fltarr(3)

  restore, 'reionfiletemplate.sav'
  reiondata = read_ascii('set36/reion.out', template=reionfiletemplate)
  redshift = reiondata.z
  gpi = reiondata.field04

  readcol, '../data/gammapi_cafg.dat', x, y, dy, /silent 
  x[0] = 0.0 
  oplot, x, y, psym=8, color=2
  oploterror, x, y, dy, errcolor=2, psym=3

  legend, ['Faucher-Giguere 08'], linestyle=[0], color=[2], psym=[8], /bottom, charsize=1

  device, /close_file
  set_plot, 'X'

END
