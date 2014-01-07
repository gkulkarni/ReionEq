
PRO abr, row, z

  ylbound = 1.0e-6
  yubound = 1.0e0

  set_plot, 'ps'
  device, filename='dla_dist.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

  ;; window, xsize=1000, ysize=1000
  !P.multi = [0, 1, 2, 0, 1]
  ;; Device, decomposed=0
  TvLCT, 255, 0, 0, 2 
  TvLCT, 0, 127, 255, 3
  !P.charsize = 1.5
  !P.thick = 1.0
  !P.charthick = 1

  redshift_label = 'z='+strtrim(uint(z),2)

  smallh = 0.71 
  lcolumn = 86
  yrbys = 3.154e7
  cmbympc = 3.241e-25 
  hubble_0 = 1.023e-10*0.71     ; yr^-1 
  omega_nr = 0.2646
  omega_lambda = 0.734
  omega_k = 1.0-(omega_lambda+omega_nr)  
  pi = 3.14159265358979
  speed_of_light = 2.998e10     ; cm/s
  mproton = 1.672622e-27        ; mass of proton ; kg
  restore, 'halo_template.sav'

  cutoff = 1.0e-2
  
;---------------------------

@cbyfe_f 
  for i = 0, size(cs3_cbyfe2, /n_elements)-1 do begin 
     if (ns3_cbyfe[i] lt cutoff) then ns3_cbyfe[i] = -1.0 
     if (ns1_cbyfe[i] lt cutoff) then ns1_cbyfe[i] = -1.0 
  endfor

  plot, cs3_cbyfe2, ns3_cbyfe, xtitle='!6[C/Fe]', ytitle='d!E2!NN/dXd[C/Fe]', xrange=[-0.6,0.6], $
        xstyle=1, yrange=[0.0, 2.0], ystyle=1
  oplot, cs1_cbyfe2, ns1_cbyfe, color=2

@cbyfe_m 
  for i = 0, size(cs3_cbyfe2, /n_elements)-1 do begin 
     if (ns3_cbyfe[i] lt cutoff) then ns3_cbyfe[i] = -1.0 
     if (ns1_cbyfe[i] lt cutoff) then ns1_cbyfe[i] = -1.0 
  endfor

  oplot, cs3_cbyfe2, ns3_cbyfe, thick=3
  oplot, cs1_cbyfe2, ns1_cbyfe, color=2, thick=3

@cbyfe_o 
  for i = 0, size(cs3_cbyfe2, /n_elements)-1 do begin 
     ; print, i, cs3_cbyfe2[i], ns3_cbyfe[i]
     if (ns3_cbyfe[i] lt cutoff) then ns3_cbyfe[i] = -1.0 
     if (ns1_cbyfe[i] lt cutoff) then ns1_cbyfe[i] = -1.0 
  endfor

  oplot, cs3_cbyfe2, ns3_cbyfe, thick=6
  oplot, cs1_cbyfe2, ns1_cbyfe, color=2, thick=6

;---------------------------

@obysi_f 
  for i = 0, size(cs3_obysi2, /n_elements)-1 do begin 
     if (ns3_obysi[i] lt cutoff) then ns3_obysi[i] = -1.0 
     if (ns1_obysi[i] lt cutoff) then ns1_obysi[i] = -1.0 
  endfor
  plot, cs3_obysi2, ns3_obysi, xtitle='[O/Si]', ytitle='d!E2!NN/dXd[O/Si]', xrange=[-0.6,0.6], $
        xstyle=1, yrange=[0.0, 2.0], ystyle=1
  oplot, cs1_obysi2, ns1_obysi, color=2

@obysi_m 
  for i = 0, size(cs3_obysi2, /n_elements)-1 do begin 
     if (ns3_obysi[i] lt cutoff) then ns3_obysi[i] = -1.0 
     if (ns1_obysi[i] lt cutoff) then ns1_obysi[i] = -1.0 
  endfor

  oplot, cs3_obysi2, ns3_obysi, thick=3
  oplot, cs1_obysi2, ns1_obysi, color=2, thick=3

@obysi_o
  for i = 0, size(cs3_obysi2, /n_elements)-1 do begin 
     if (ns3_obysi[i] lt cutoff) then ns3_obysi[i] = -1.0 
     if (ns1_obysi[i] lt cutoff) then ns1_obysi[i] = -1.0 
  endfor
  ns3_obysi[0] = -1.0
  ns3_obysi = ts_smooth(ns3_obysi, 40) 
  oplot, cs3_obysi2, ns3_obysi, thick=6
  ns1_obysi = ts_smooth(ns1_obysi, 40) 
  oplot, cs1_obysi2, ns1_obysi, color=2, thick=6

  legend, ['t!Ddelay!N=0', 't!Ddelay!N=10!E11!Nyr', $
           't!Ddelay!N=10!E12!Nyr'], $
          linestyle=[0,0,0], thick=[1,3,6], color=[-1,-1,-1], charsize=1.1

device, /close_file
set_plot, 'X'

END

