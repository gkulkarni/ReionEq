pro plot_reion 

  openr, 1, '../reion.out'
  z = 0. 
  arg = fltarr(12) 
  rs = fltarr(100)
  rs_log = fltarr(100)
  q = fltarr(100) 
  gpi = fltarr(100) 
  i = 0 
  while not eof(1) do begin 
     readf, 1, z, arg 
     rs[i] = z 
     ;; rs_log[i] = alog10(z)
     q[i] = arg[0]
     ;; gpi[i] = alog10(arg[2]/1.0e-10)
     gpi[i] = arg[2]/1.0e-10
     i = i + 1 
  endwhile
  close, 1

  ; HII region filling factor.
  set_plot, 'ps'
  device, filename='q.ps'
  !p.font = -1
  plot, rs, q, xtitle = '!6redshift', ytitle = 'filling factor', xrange=[0.,30.]
  device, /close 

  device, filename='gpi.ps'
  !p.font = -1
  plot, rs, gpi, xtitle = '!6redshift', ytitle = 'photoionization rate', xrange=[3.,20.], yrange=[0.01,5.], /xlog, /ylog
  device, /close 
  
end
