
PRO sfr, set

  window, xsize=1000, ysize=1000
  Device, decomposed=0
  TvLCT, 255, 0, 0, 2 
  !P.charsize = 2

  readcol, '../data/sfrsfr.data', rs, rserr, sf, sferr, /silent  
  sfe = 10.0^sf 
  plotsym, 0, 0.5, /FILL
  plot, rs, sfe, psym=8, /ylog, /xlog, xrange=[5.0e-1,100], $
        ytitle='SFR (M!D!9n!N!X yr!E-1!N Mpc!E-3 !N)', xstyle=1, yrange=[1.0e-4,1], $
        xtitle='!6z', /nodata 
  oplot, rs, sfe, psym=8, color=2
  yup = 10.0^(sf+sferr*0.5)
  ylo = 10.0^(sf-sferr*0.5)
  dy = yup-ylo
  plot_err, rs, sfe, dy, dx1=rserr, color=2

  restore, 'sfrfiletemplate.sav'
  sfrdata = read_ascii(set + '/sfr.out', template=sfrfiletemplate)
  sfr_tot = sfrdata.total_sfr
  z = sfrdata.redshift
  oplot, z, sfr_tot 
  sfr_pop2 = sfrdata.pop2_sfr
  sfr_pop3 = sfrdata.pop3_sfr
  n = size(sfr_pop3, /n_elements)
  oplot, z, sfr_pop3, linestyle=5

END

