
; File: abr2.pro
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; This is an important piece of code. It creates the two central
; result figures of our paper.  It plots various abundance ratios by
; reading the results of reion.  It makes two kinds of plot for each
; abundance ratio: (1) ratio vs. halo mass, and (2) the DLA
; distribution function vs. abundance ratio.  It takes the redshift as
; a command line argument. 

PRO abr2, row, z

  ylbound = 1.0e-6
  yubound = 1.0e0

;; set_plot, 'ps'
;; device, filename='abr2_z6.ps', xsize=7.0, ysize=10.0, /inches, color=1, yoffset=1.0

  window, xsize=1000, ysize=1000
  !P.multi = [0, 2, 4, 0, 1]
  Device, decomposed=0
  TvLCT, 255, 0, 0, 2 
  TvLCT, 0, 127, 255, 3
  !P.charsize = 1.5
  !P.thick = 1.0
  !y.omargin = [2,6]

  redshift_label = 'z='+strtrim(uint(z),2)

;----------------------------------------

  smallh = 0.71 
  lcolumn = 95
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

;----------------------------------------

  c_data2 = read_ascii('set95/halos_c.out', template=stars_template)
  fe_data2 = read_ascii('set95/halos_fe.out', template=stars_template)
  m_data2 = read_ascii('set95/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set95/halos_gas.out', template=stars_template)

  c2 = c_data2.field001[lcolumn:*,row]
  fe2 = fe_data2.field001[lcolumn:*,row]
  m2 = m_data2.field001[lcolumn:*,row]
  n2 = n_data2.field001[lcolumn:*,row]
  gas2 = gas_data2.field001[lcolumn:*,row]

  array_size = size(c2)
  array_length = array_size[1] 

  cbyfe2 = fltarr(array_length)
  febyh2 = fltarr(array_length)
  cross_section = fltarr(array_length)
  amplitude = fltarr(array_length)
  dndx = fltarr(array_length)
  foo = fltarr(array_length) 

  sigma0 = 40.0                 ; kpc^2
  m0 = 10.0^(-0.5)              ; 10^10 M_solar 
  alpha = 0.2                   ; dimensionless 
  
  ;; sigma0 = 48.9 
  ;; m0 = 10.0^(-0.5)
  ;; alpha = 1.0

  ;; sigma0 = 40.0
  ;; m0 = 1.0
  ;; alpha = 0.5

  for i = 0, array_length-1 do begin 

     cnst = -0.05 

     if (i eq 0) then begin 
        dm = 0.01 * m2[i]
     endif else begin 
        dm = m2[i]-m2[i-1]
     endelse

     mH2 = gas2[i]*0.71
     
     if (fe2[i] eq 0.0) then begin 
        febyh2[i] = -10.0
     endif else begin 
        febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
     endelse
     
     if (c2[i] eq 0.0) then begin 
        cbyfe2[i] = -10.0
     endif else begin 
        cbyfe2[i] = alog10(abs(c2[i]/fe2[i])) - 0.41
     endelse

                                ; We need to find mass of a halo at
                                ; z=3 that has the same circular
                                ; velocity as this halo. 
     omega_mz = omega_nr*(1.0+z)^3/(omega_nr*(1.0+z)^3+omega_lambda+omega_k*(1.0+z)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     vc = 23.4 * (m2[i]*smallh/1.0e-2)^(1.0/3.0) $
          * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) $
          * ((1.0+z)/10.0)^0.5  ; km/s 

     rvir = 0.784 * (m2[i]/(1.0e-2/smallh))^(1.0/3.0) *$
            ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(-1.0/3.0) *$
            (10.0/(1.0+z)) / smallh ; kpc 

     omega_mz = omega_nr*(4.0)^3/(omega_nr*(4.0)^3+omega_lambda+omega_k*(4.0)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     m_z3 = vc/(((4.0)/10.0)^0.5 * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) * (smallh/1.0e-2)^(1.0/3.0) * 23.4)
     m_z3 = m_z3^3 

                                ; cross_section[i] = sigma0 * (m2[i]/m0)^2 * (1.0 + m2[i]/m0)^(alpha-2.0) ; kpc^2 
     cross_section[i] = sigma0 * (m_z3/m0)^2 * (1.0 + m_z3/m0)^(alpha-2.0) ; kpc^2 
     amplitude[i] = cross_section[i]*n2[i]*1.0e-6                          ; dimensionless (1.0e-6 is kpc^2/Mpc^2)

                                ; dndx is a misnomer here! What I mean
                                ; is (d^2n/dXdp . dp).  See discussion
                                ; around Eqn. 3 of Pontzen et al. (2008)
                                ; dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/(hubble_0*(1.0+z)^3) ; dimensionless 
     dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/hubble_0 ; dimensionless 
     foo[i] = dndx[i] 

  endfor

  

  deepee = fltarr(array_length) 

  mean = 0.0 
  ldla = 0.0 
  for i = 0, array_length-1 do begin 
     if (i eq 0) then begin 
        dp = abs(cbyfe2[i])
     endif else begin 
        dp = abs(cbyfe2[i] - cbyfe2[i-1])
     endelse
     deepee[i] = dp 
     if (dp EQ 0.0) then dp = 0.001
     dndx[i] = dndx[i]/dp
     ldla = ldla + (dndx[i] * dp) 
     mean = mean + (dndx[i] * cbyfe2[i] * dp)
  endfor 

  a = float(ylbound)
  b = float(yubound) 

  m2 = m2*1.0e10 
; plot, m2, cbyfe2, /xlog, xrange=[1.0e8, 1.0e12], psym=-6, yrange=[-0.2,1.0]
  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 
  norm = 0.0 

  for i = 0, array_length-1 do begin 
     norm = norm + foo[i] 
     dndx[i] = foo[i] / dp3[i] 
  endfor 

; print, 'norm=', norm 

  m3 = m2 
  cbyfe3 = cbyfe2 
  dndx3 = dndx

; plot, cbyfe2, dndx, psym=-6

  mean = mean/ldla

  print, 'ldla=', ldla 
; print, z, ldla
; return 
  print, 'mean [c/fe] chosen set=', mean

;  return 

;--------------------------------------

  c_data2 = read_ascii('set49/halos_c.out', template=stars_template)
  fe_data2 = read_ascii('set49/halos_fe.out', template=stars_template)
  m_data2 = read_ascii('set49/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set49/halos_gas.out', template=stars_template)

  c2 = c_data2.field001[lcolumn:*,row]
  fe2 = fe_data2.field001[lcolumn:*,row]
  m2 = m_data2.field001[lcolumn:*,row]
  n2 = n_data2.field001[lcolumn:*,row]
  gas2 = gas_data2.field001[lcolumn:*,row]

  array_size = size(c2)
  array_length = array_size[1] 

  cbyfe2 = fltarr(array_length)
  febyh2 = fltarr(array_length)
  cross_section = fltarr(array_length)
  amplitude = fltarr(array_length)
  dndx = fltarr(array_length)
  foo = fltarr(array_length) 

;  sigma0 = 40.0                 ; kpc^2
;  m0 = 10.0^(-0.5)              ; 10^10 M_solar 
;  alpha = 0.2                   ; dimensionless 

  for i = 0, array_length-1 do begin 

     cnst = -0.05 

     if (i eq 0) then begin 
        dm = 0.01 * m2[i]
     endif else begin 
        dm = m2[i]-m2[i-1]
     endelse

     mH2 = gas2[i]*0.71
     
     if (fe2[i] eq 0.0) then begin 
        febyh2[i] = -10.0
     endif else begin 
        febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
     endelse
;     print, i, febyh2[i] 

     
     if (c2[i] eq 0.0) then begin 
        cbyfe2[i] = -10.0
     endif else begin 
        cbyfe2[i] = alog10(abs(c2[i]/fe2[i])) - 0.41
     endelse

                                ; We need to find mass of a halo at
                                ; z=3 that has the same circular
                                ; velocity as this halo. 
     omega_mz = omega_nr*(1.0+z)^3/(omega_nr*(1.0+z)^3+omega_lambda+omega_k*(1.0+z)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     vc = 23.4 * (m2[i]*smallh/1.0e-2)^(1.0/3.0) $
          * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) $
          * ((1.0+z)/10.0)^0.5  ; km/s 

     omega_mz = omega_nr*(4.0)^3/(omega_nr*(4.0)^3+omega_lambda+omega_k*(4.0)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     m_z3 = vc/(((4.0)/10.0)^0.5 * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) * (smallh/1.0e-2)^(1.0/3.0) * 23.4)
     m_z3 = m_z3^3 

                                ; cross_section[i] = sigma0 * (m2[i]/m0)^2 * (1.0 + m2[i]/m0)^(alpha-2.0) ; kpc^2 
     cross_section[i] = sigma0 * (m_z3/m0)^2 * (1.0 + m_z3/m0)^(alpha-2.0) ; kpc^2 
     amplitude[i] = cross_section[i]*n2[i]*1.0e-6                          ; dimensionless (1.0e-6 is kpc^2/Mpc^2)

                                ; dndx is a misnomer here! What I mean
                                ; is (d^2n/dXdp . dp).  See discussion
                                ; around Eqn. 3 of Pontzen et al. (2008)
                                ; dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/(hubble_0*(1.0+z)^3) ; dimensionless 
     dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/hubble_0 ; dimensionless 
     foo[i] = dndx[i] 

  endfor

  deepee = fltarr(array_length) 

  mean = 0.0 
  ldla = 0.0 
  for i = 0, array_length-1 do begin 

     if (i eq 0) then begin 
        dp = abs(cbyfe2[i])
     endif else begin 
        dp = abs(cbyfe2[i] - cbyfe2[i-1])
     endelse

     deepee[i] = dp 

     if (dp EQ 0.0) then dp = 0.001

     dndx[i] = dndx[i]/dp

     ldla = ldla + (dndx[i] * dp) 

     mean = mean + (dndx[i] * cbyfe2[i] * dp)

  endfor 

  a = float(ylbound)
  b = float(yubound) 

  m2 = m2*1.0e10 
; plot, m2, cbyfe2, /xlog, xrange=[1.0e8, 1.0e12], psym=-6, yrange=[-0.2,1.0]
  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 
  norm = 0.0 

  for i = 0, array_length-1 do begin 

     norm = norm + foo[i] 

     dndx[i] = foo[i] / dp3[i] 

  endfor 

  print, 'norm=', norm 

  m1 = m2 
  cbyfe1 = cbyfe2 
  dndx1 = dndx


; plot, cbyfe2, dndx, psym=-6

  mean = mean/ldla

  print, 'ldla=', ldla 
  print, 'mean[c/fe]=', mean


;----------------------------------------

  c_data2 = read_ascii('set50/halos_c.out', template=stars_template)
  fe_data2 = read_ascii('set50/halos_fe.out', template=stars_template)
  m_data2 = read_ascii('set50/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set50/halos_gas.out', template=stars_template)

  c2 = c_data2.field001[lcolumn:*,row]
  fe2 = fe_data2.field001[lcolumn:*,row]
  m2 = m_data2.field001[lcolumn:*,row]
  n2 = n_data2.field001[lcolumn:*,row]
  gas2 = gas_data2.field001[lcolumn:*,row]

  array_size = size(c2)
  array_length = array_size[1] 

  cbyfe2 = fltarr(array_length)
  febyh2 = fltarr(array_length)
  cross_section = fltarr(array_length)
  amplitude = fltarr(array_length)
  dndx = fltarr(array_length)
  foo = fltarr(array_length) 

;  sigma0 = 40.0                 ; kpc^2
;  m0 = 10.0^(-0.5)              ; 10^10 M_solar 
;  alpha = 0.2                   ; dimensionless 

  for i = 0, array_length-1 do begin 

     cnst = -0.05 

     if (i eq 0) then begin 
        dm = 0.01 * m2[i]
     endif else begin 
        dm = m2[i]-m2[i-1]
     endelse

     mH2 = gas2[i]*0.71
     
     if (fe2[i] eq 0.0) then begin 
        febyh2[i] = -10.0
     endif else begin 
        febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
     endelse
     
     if (c2[i] eq 0.0) then begin 
        cbyfe2[i] = -10.0
     endif else begin 
        cbyfe2[i] = alog10(abs(c2[i]/fe2[i])) - 0.41
     endelse

                                ; We need to find mass of a halo at
                                ; z=3 that has the same circular
                                ; velocity as this halo. 
     omega_mz = omega_nr*(1.0+z)^3/(omega_nr*(1.0+z)^3+omega_lambda+omega_k*(1.0+z)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     vc = 23.4 * (m2[i]*smallh/1.0e-2)^(1.0/3.0) $
          * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) $
          * ((1.0+z)/10.0)^0.5  ; km/s 

     omega_mz = omega_nr*(4.0)^3/(omega_nr*(4.0)^3+omega_lambda+omega_k*(4.0)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     m_z3 = vc/(((4.0)/10.0)^0.5 * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) * (smallh/1.0e-2)^(1.0/3.0) * 23.4)
     m_z3 = m_z3^3 

                                ; cross_section[i] = sigma0 * (m2[i]/m0)^2 * (1.0 + m2[i]/m0)^(alpha-2.0) ; kpc^2 
     cross_section[i] = sigma0 * (m_z3/m0)^2 * (1.0 + m_z3/m0)^(alpha-2.0) ; kpc^2 
     amplitude[i] = cross_section[i]*n2[i]*1.0e-6                          ; dimensionless (1.0e-6 is kpc^2/Mpc^2)

                                ; dndx is a misnomer here! What I mean
                                ; is (d^2n/dXdp . dp).  See discussion
                                ; around Eqn. 3 of Pontzen et al. (2008)
                                ; dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/(hubble_0*(1.0+z)^3) ; dimensionless 
     dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/hubble_0 ; dimensionless 
     foo[i] = dndx[i] 

  endfor

  deepee = fltarr(array_length) 

  mean = 0.0 
  ldla = 0.0 
  for i = 0, array_length-1 do begin 

     if (i eq 0) then begin 
        dp = abs(cbyfe2[i])
     endif else begin 
        dp = abs(cbyfe2[i] - cbyfe2[i-1])
     endelse

     deepee[i] = dp 

     if (dp EQ 0.0) then dp = 0.001

     dndx[i] = dndx[i]/dp

     ldla = ldla + (dndx[i] * dp) 

     mean = mean + (dndx[i] * cbyfe2[i] * dp)

  endfor 

  a = float(ylbound)
  b = float(yubound) 

  m2 = m2*1.0e10 

  ms = rebin(m2, array_length/2)
  cs = rebin(cbyfe2, array_length/2)
  ms1 = rebin(m1, array_length/2)
  cs1 = rebin(cbyfe1, array_length/2)
  ms3 = rebin(m3, array_length/2)
  cs3 = rebin(cbyfe3, array_length/2)

;; plot, ms3, cs3, /xlog, xrange=[1.0e8, 1.0e12], psym=-6, yrange=[-0.2,1.0], xtitle='halo mass', ytitle='[C/Fe]', symsize=0.5
;; oplot, ms1, cs1, color=2, psym=-6, symsize=0.5
;; oplot, ms, cs, color=3, psym=-6, symsize=0.5

  ms_cbyfe = ms
  cs_cbyfe = cs
  ms1_cbyfe = ms1
  cs1_cbyfe = cs1
  ms3_cbyfe = ms3
  cs3_cbyfe = cs3 

  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 

  norm = 0.0 
  for i = 0, array_length-1 do begin 
     norm = norm + foo[i] 
     dndx[i] = foo[i] / dp3[i] 
  endfor 
  print, 'norm=', norm 

  ns = rebin(dndx, array_length/2) 
  cs = rebin(cbyfe2, array_length/2)
  ns1 = rebin(dndx1, array_length/2) 
  cs1 = rebin(cbyfe1, array_length/2)
  ns3 = rebin(dndx3, array_length/2) 
  cs3 = rebin(cbyfe3, array_length/2)

;; plot, cs3, ns3, psym=-6, xtitle='[C/Fe]', ytitle='d!E2!NN/dXd[C/Fe]', symsize=0.5, xrange=[-0.2,1.0]
;; oplot, cs1, ns1, psym=-6, color=2, symsize=0.5 
;; oplot, cs, ns, psym=6, color=3, symsize=0.5

  ns_cbyfe = ns
  cs_cbyfe2 = cs
  ns1_cbyfe = ns1
  cs1_cbyfe2 = cs1
  ns3_cbyfe = ns3
  cs3_cbyfe2 = cs3

  mean = mean/ldla
  print, 'ldla=', ldla 
  print, 'mean[c/fe]=', mean

;; >>>>>1<<<<<

;-----------------

; IMPORTANT below:
; C = O 
; Fe = Si 

  c_data2 = read_ascii('set95/halos_o.out', template=stars_template)
  fe_data2 = read_ascii('set95/halos_Si.out', template=stars_template)
  m_data2 = read_ascii('set95/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set95/halos_gas.out', template=stars_template)

  c2 = c_data2.field001[lcolumn:*,row]
  fe2 = fe_data2.field001[lcolumn:*,row]
  m2 = m_data2.field001[lcolumn:*,row]
  n2 = n_data2.field001[lcolumn:*,row]
  gas2 = gas_data2.field001[lcolumn:*,row]

  array_size = size(c2)
  array_length = array_size[1] 

  cbyfe2 = fltarr(array_length)
  febyh2 = fltarr(array_length)
  cross_section = fltarr(array_length)
  amplitude = fltarr(array_length)
  dndx = fltarr(array_length)
  foo = fltarr(array_length) 

;  sigma0 = 40.0                 ; kpc^2
;  m0 = 10.0^(-0.5)              ; 10^10 M_solar 
;  alpha = 0.2                   ; dimensionless 

  for i = 0, array_length-1 do begin 

     cnst = -0.05 

     if (i eq 0) then begin 
        dm = 0.01 * m2[i]
     endif else begin 
        dm = m2[i]-m2[i-1]
     endelse

     mH2 = gas2[i]*0.71
     
     if (fe2[i] eq 0.0) then begin 
        febyh2[i] = -10.0
     endif else begin 
        febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
     endelse
     
     if (c2[i] eq 0.0) then begin 
        cbyfe2[i] = -10.0
     endif else begin 
        cbyfe2[i] = alog10(abs(c2[i]/fe2[i])) - 1.17
     endelse

                                ; We need to find mass of a halo at
                                ; z=3 that has the same circular
                                ; velocity as this halo. 
     omega_mz = omega_nr*(1.0+z)^3/(omega_nr*(1.0+z)^3+omega_lambda+omega_k*(1.0+z)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     vc = 23.4 * (m2[i]*smallh/1.0e-2)^(1.0/3.0) $
          * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) $
          * ((1.0+z)/10.0)^0.5  ; km/s 

     omega_mz = omega_nr*(4.0)^3/(omega_nr*(4.0)^3+omega_lambda+omega_k*(4.0)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     m_z3 = vc/(((4.0)/10.0)^0.5 * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) * (smallh/1.0e-2)^(1.0/3.0) * 23.4)
     m_z3 = m_z3^3 

     rvir = 0.784 * (m2[i]/(1.0e-2/smallh))^(1.0/3.0) *$
            ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(-1.0/3.0) *$
            (10.0/(1.0+z)) / smallh ; kpc 
;     print, i, m2[i], 3.0*double(c2[i])*1.0e10*1.98e30/(2.0*pi*16.0*double(mproton)*0.8*0.8*(double(rvir)*3.09e21)^2)


                                ; cross_section[i] = sigma0 * (m2[i]/m0)^2 * (1.0 + m2[i]/m0)^(alpha-2.0) ; kpc^2 
     cross_section[i] = sigma0 * (m_z3/m0)^2 * (1.0 + m_z3/m0)^(alpha-2.0) ; kpc^2 
     amplitude[i] = cross_section[i]*n2[i]*1.0e-6                          ; dimensionless (1.0e-6 is kpc^2/Mpc^2)

                                ; dndx is a misnomer here! What I mean
                                ; is (d^2n/dXdp . dp).  See discussion
                                ; around Eqn. 3 of Pontzen et al. (2008)
                                ; dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/(hubble_0*(1.0+z)^3) ; dimensionless 
     dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/hubble_0 ; dimensionless 
     foo[i] = dndx[i] 

  endfor

  deepee = fltarr(array_length) 

  mean = 0.0 
  ldla = 0.0 
  for i = 0, array_length-1 do begin 

     if (i eq 0) then begin 
        dp = abs(cbyfe2[i])
     endif else begin 
        dp = abs(cbyfe2[i] - cbyfe2[i-1])
     endelse

     deepee[i] = dp 

     if (dp EQ 0.0) then dp = 0.001

     dndx[i] = dndx[i]/dp

     ldla = ldla + (dndx[i] * dp) 

     mean = mean + (dndx[i] * cbyfe2[i] * dp)

  endfor 

  a = float(ylbound)
  b = float(yubound) 

  m2 = m2*1.0e10 
; plot, m2, cbyfe2, /xlog, xrange=[1.0e8, 1.0e12], psym=-6, yrange=[-0.2,1.0]
  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 
  norm = 0.0 

  for i = 0, array_length-1 do begin 

     norm = norm + foo[i] 

     dndx[i] = foo[i] / dp3[i] 

  endfor 

  print, 'norm=', norm 

  m3 = m2 
  cbyfe3 = cbyfe2 
  dndx3 = dndx


; plot, cbyfe2, dndx, psym=-6

  mean = mean/ldla

  print, 'ldla=', ldla 
  print, 'mean[o/si] chosen set=', mean

;--------------------------------------

  c_data2 = read_ascii('set49/halos_o.out', template=stars_template)
  fe_data2 = read_ascii('set49/halos_Si.out', template=stars_template)
  m_data2 = read_ascii('set49/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set49/halos_gas.out', template=stars_template)

  c2 = c_data2.field001[lcolumn:*,row]
  fe2 = fe_data2.field001[lcolumn:*,row]
  m2 = m_data2.field001[lcolumn:*,row]
  n2 = n_data2.field001[lcolumn:*,row]
  gas2 = gas_data2.field001[lcolumn:*,row]

  array_size = size(c2)
  array_length = array_size[1] 

  cbyfe2 = fltarr(array_length)
  febyh2 = fltarr(array_length)
  cross_section = fltarr(array_length)
  amplitude = fltarr(array_length)
  dndx = fltarr(array_length)
  foo = fltarr(array_length) 

;  sigma0 = 40.0                 ; kpc^2
;  m0 = 10.0^(-0.5)              ; 10^10 M_solar 
;  alpha = 0.2                   ; dimensionless 

  for i = 0, array_length-1 do begin 

     cnst = -0.05 

     if (i eq 0) then begin 
        dm = 0.01 * m2[i]
     endif else begin 
        dm = m2[i]-m2[i-1]
     endelse

     mH2 = gas2[i]*0.71
     
     if (fe2[i] eq 0.0) then begin 
        febyh2[i] = -10.0
     endif else begin 
        febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
     endelse
     
     if (c2[i] eq 0.0) then begin 
        cbyfe2[i] = -10.0
     endif else begin 
        cbyfe2[i] = alog10(abs(c2[i]/fe2[i])) - 1.17
     endelse

                                ; We need to find mass of a halo at
                                ; z=3 that has the same circular
                                ; velocity as this halo. 
     omega_mz = omega_nr*(1.0+z)^3/(omega_nr*(1.0+z)^3+omega_lambda+omega_k*(1.0+z)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     vc = 23.4 * (m2[i]*smallh/1.0e-2)^(1.0/3.0) $
          * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) $
          * ((1.0+z)/10.0)^0.5  ; km/s 

     omega_mz = omega_nr*(4.0)^3/(omega_nr*(4.0)^3+omega_lambda+omega_k*(4.0)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     m_z3 = vc/(((4.0)/10.0)^0.5 * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) * (smallh/1.0e-2)^(1.0/3.0) * 23.4)
     m_z3 = m_z3^3 

                                ; cross_section[i] = sigma0 * (m2[i]/m0)^2 * (1.0 + m2[i]/m0)^(alpha-2.0) ; kpc^2 
     cross_section[i] = sigma0 * (m_z3/m0)^2 * (1.0 + m_z3/m0)^(alpha-2.0) ; kpc^2 
     amplitude[i] = cross_section[i]*n2[i]*1.0e-6                          ; dimensionless (1.0e-6 is kpc^2/Mpc^2)

                                ; dndx is a misnomer here! What I mean
                                ; is (d^2n/dXdp . dp).  See discussion
                                ; around Eqn. 3 of Pontzen et al. (2008)
                                ; dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/(hubble_0*(1.0+z)^3) ; dimensionless 
     dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/hubble_0 ; dimensionless 
     foo[i] = dndx[i] 

  endfor

  deepee = fltarr(array_length) 

  mean = 0.0 
  ldla = 0.0 
  for i = 0, array_length-1 do begin 

     if (i eq 0) then begin 
        dp = abs(cbyfe2[i])
     endif else begin 
        dp = abs(cbyfe2[i] - cbyfe2[i-1])
     endelse

     deepee[i] = dp 

     if (dp EQ 0.0) then dp = 0.001

     dndx[i] = dndx[i]/dp

     ldla = ldla + (dndx[i] * dp) 

     mean = mean + (dndx[i] * cbyfe2[i] * dp)

  endfor 

  a = float(ylbound)
  b = float(yubound) 

  m2 = m2*1.0e10 
; plot, m2, cbyfe2, /xlog, xrange=[1.0e8, 1.0e12], psym=-6, yrange=[-0.2,1.0]
  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 
  norm = 0.0 

  for i = 0, array_length-1 do begin 

     norm = norm + foo[i] 

     dndx[i] = foo[i] / dp3[i] 

  endfor 

  print, 'norm=', norm 

  m1 = m2 
  cbyfe1 = cbyfe2 
  dndx1 = dndx


; plot, cbyfe2, dndx, psym=-6

  mean = mean/ldla

  print, 'ldla=', ldla 
  print, 'mean[o/si]=', mean


;----------------------------------------

  c_data2 = read_ascii('set50/halos_o.out', template=stars_template)
  fe_data2 = read_ascii('set50/halos_Si.out', template=stars_template)
  m_data2 = read_ascii('set50/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set50/halos_gas.out', template=stars_template)

  c2 = c_data2.field001[lcolumn:*,row]
  fe2 = fe_data2.field001[lcolumn:*,row]
  m2 = m_data2.field001[lcolumn:*,row]
  n2 = n_data2.field001[lcolumn:*,row]
  gas2 = gas_data2.field001[lcolumn:*,row]

  array_size = size(c2)
  array_length = array_size[1] 

  cbyfe2 = fltarr(array_length)
  febyh2 = fltarr(array_length)
  cross_section = fltarr(array_length)
  amplitude = fltarr(array_length)
  dndx = fltarr(array_length)
  foo = fltarr(array_length) 

;  sigma0 = 40.0                 ; kpc^2
;  m0 = 10.0^(-0.5)              ; 10^10 M_solar 
;  alpha = 0.2                   ; dimensionless 

  for i = 0, array_length-1 do begin 

     cnst = -0.05 

     if (i eq 0) then begin 
        dm = 0.01 * m2[i]
     endif else begin 
        dm = m2[i]-m2[i-1]
     endelse

     mH2 = gas2[i]*0.71
     
     if (fe2[i] eq 0.0) then begin 
        febyh2[i] = -10.0
     endif else begin 
        febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
     endelse
     
     if (c2[i] eq 0.0) then begin 
        cbyfe2[i] = -10.0
     endif else begin 
        cbyfe2[i] = alog10(abs(c2[i]/fe2[i])) - 1.17
     endelse

                                ; We need to find mass of a halo at
                                ; z=3 that has the same circular
                                ; velocity as this halo. 
     omega_mz = omega_nr*(1.0+z)^3/(omega_nr*(1.0+z)^3+omega_lambda+omega_k*(1.0+z)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     vc = 23.4 * (m2[i]*smallh/1.0e-2)^(1.0/3.0) $
          * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) $
          * ((1.0+z)/10.0)^0.5  ; km/s 

     omega_mz = omega_nr*(4.0)^3/(omega_nr*(4.0)^3+omega_lambda+omega_k*(4.0)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     m_z3 = vc/(((4.0)/10.0)^0.5 * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) * (smallh/1.0e-2)^(1.0/3.0) * 23.4)
     m_z3 = m_z3^3 

                                ; cross_section[i] = sigma0 * (m2[i]/m0)^2 * (1.0 + m2[i]/m0)^(alpha-2.0) ; kpc^2 
     cross_section[i] = sigma0 * (m_z3/m0)^2 * (1.0 + m_z3/m0)^(alpha-2.0) ; kpc^2 
     amplitude[i] = cross_section[i]*n2[i]*1.0e-6                          ; dimensionless (1.0e-6 is kpc^2/Mpc^2)

                                ; dndx is a misnomer here! What I mean
                                ; is (d^2n/dXdp . dp).  See discussion
                                ; around Eqn. 3 of Pontzen et al. (2008)
                                ; dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/(hubble_0*(1.0+z)^3) ; dimensionless 
     dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/hubble_0 ; dimensionless 
     foo[i] = dndx[i] 

  endfor

  deepee = fltarr(array_length) 

  mean = 0.0 
  ldla = 0.0 
  for i = 0, array_length-1 do begin 

     if (i eq 0) then begin 
        dp = abs(cbyfe2[i])
     endif else begin 
        dp = abs(cbyfe2[i] - cbyfe2[i-1])
     endelse

     deepee[i] = dp 

     if (dp EQ 0.0) then dp = 0.001

     dndx[i] = dndx[i]/dp

     ldla = ldla + (dndx[i] * dp) 

     mean = mean + (dndx[i] * cbyfe2[i] * dp)

  endfor 

  a = float(ylbound)
  b = float(yubound) 

  m2 = m2*1.0e10 

  ms = rebin(m2, array_length/2)
  cs = rebin(cbyfe2, array_length/2)
  ms1 = rebin(m1, array_length/2)
  cs1 = rebin(cbyfe1, array_length/2)
  ms3 = rebin(m3, array_length/2)
  cs3 = rebin(cbyfe3, array_length/2)

;; plot, ms3, cs3, /xlog, xrange=[1.0e8, 1.0e12], psym=-6, yrange=[-0.5,1.0], xtitle='halo mass', ytitle='[O/Si]', symsize=0.5
;; oplot, ms1, cs1, color=2, psym=-6, symsize=0.5
;; oplot, ms, cs, color=3, psym=-6, symsize=0.5

  ms_obysi = ms
  cs_obysi = cs
  ms1_obysi = ms1
  cs1_obysi = cs1
  ms3_obysi = ms3
  cs3_obysi = cs3

  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 

  norm = 0.0 
  for i = 0, array_length-1 do begin 
     norm = norm + foo[i] 
     dndx[i] = foo[i] / dp3[i] 
  endfor 
  print, 'norm=', norm 

  n = 50 
  t = fltarr(n)
  dx = (max(cbyfe2)-min(cbyfe2))/n 
  xmin = min(cbyfe2) 
  for i = 0, n-1 do begin 
     t[i] = xmin + dx*i 
  endfor 

  ns = rebin(dndx, array_length/2) 
  cs = rebin(cbyfe2, array_length/2)

;plot, cbyfe3, dndx3, psym=-6, xtitle='[O/Si]', ytitle='d!E2!NN/dXd[O/Si]', xrange=[-0.5,1.0], symsize=0.5

  ns1 = rebin(dndx1, array_length/2) 
  cs1 = rebin(cbyfe1, array_length/2)

;; oplot, cs1, ns1, psym=-6, color=2, symsize=0.5
;; oplot, cs, ns, psym=6, color=3, symsize=0.5

  ns_obysi = ns
  cs_obysi2 = cs
  ns1_obysi = ns1
  cs1_obysi2 = cs1
  ns3_obysi = dndx3
  cs3_obysi2 = cbyfe3

  mean = mean/ldla
  print, 'ldla=', ldla 
  print, 'mean[o/si]=', mean

;; >>>>>2<<<<<

;----------------------------------------


; IMPORTANT below:
; C = Zn 
; Fe = Fe

  c_data2 = read_ascii('set95/halos_Zn.out', template=stars_template)
  fe_data2 = read_ascii('set95/halos_fe.out', template=stars_template)
  m_data2 = read_ascii('set95/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set95/halos_gas.out', template=stars_template)

  c2 = c_data2.field001[lcolumn:*,row]
  fe2 = fe_data2.field001[lcolumn:*,row]
  m2 = m_data2.field001[lcolumn:*,row]
  n2 = n_data2.field001[lcolumn:*,row]
  gas2 = gas_data2.field001[lcolumn:*,row]

  array_size = size(c2)
  array_length = array_size[1] 

  cbyfe2 = fltarr(array_length)
  febyh2 = fltarr(array_length)
  cross_section = fltarr(array_length)
  amplitude = fltarr(array_length)
  dndx = fltarr(array_length)
  foo = fltarr(array_length) 

;  sigma0 = 40.0                 ; kpc^2
;  m0 = 10.0^(-0.5)              ; 10^10 M_solar 
;  alpha = 0.2                   ; dimensionless 

  for i = 0, array_length-1 do begin 

     cnst = -0.05 

     if (i eq 0) then begin 
        dm = 0.01 * m2[i]
     endif else begin 
        dm = m2[i]-m2[i-1]
     endelse

     mH2 = gas2[i]*0.71
     
     if (fe2[i] eq 0.0) then begin 
        febyh2[i] = -10.0
     endif else begin 
        febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
     endelse
     
     if (c2[i] eq 0.0) then begin 
        cbyfe2[i] = -10.0
     endif else begin 
        cbyfe2[i] = alog10(abs(c2[i]/fe2[i])) + 3.07
     endelse

                                ; We need to find mass of a halo at
                                ; z=3 that has the same circular
                                ; velocity as this halo. 
     omega_mz = omega_nr*(1.0+z)^3/(omega_nr*(1.0+z)^3+omega_lambda+omega_k*(1.0+z)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     vc = 23.4 * (m2[i]*smallh/1.0e-2)^(1.0/3.0) $
          * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) $
          * ((1.0+z)/10.0)^0.5  ; km/s 

     omega_mz = omega_nr*(4.0)^3/(omega_nr*(4.0)^3+omega_lambda+omega_k*(4.0)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     m_z3 = vc/(((4.0)/10.0)^0.5 * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) * (smallh/1.0e-2)^(1.0/3.0) * 23.4)
     m_z3 = m_z3^3 

                                ; cross_section[i] = sigma0 * (m2[i]/m0)^2 * (1.0 + m2[i]/m0)^(alpha-2.0) ; kpc^2 
     cross_section[i] = sigma0 * (m_z3/m0)^2 * (1.0 + m_z3/m0)^(alpha-2.0) ; kpc^2 
     amplitude[i] = cross_section[i]*n2[i]*1.0e-6                          ; dimensionless (1.0e-6 is kpc^2/Mpc^2)

                                ; dndx is a misnomer here! What I mean
                                ; is (d^2n/dXdp . dp).  See discussion
                                ; around Eqn. 3 of Pontzen et al. (2008)
                                ; dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/(hubble_0*(1.0+z)^3) ; dimensionless 
     dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/hubble_0 ; dimensionless 
     foo[i] = dndx[i] 
     
  endfor

  deepee = fltarr(array_length) 

  mean = 0.0 
  ldla = 0.0 
  for i = 0, array_length-1 do begin 

     if (i eq 0) then begin 
        dp = abs(cbyfe2[i])
     endif else begin 
        dp = abs(cbyfe2[i] - cbyfe2[i-1])
     endelse

     deepee[i] = dp 

     if (dp EQ 0.0) then dp = 0.001

     dndx[i] = dndx[i]/dp

     ldla = ldla + (dndx[i] * dp) 

     mean = mean + (dndx[i] * cbyfe2[i] * dp)

  endfor 

  a = float(ylbound)
  b = float(yubound) 

  m2 = m2*1.0e10 
; plot, m2, cbyfe2, /xlog, xrange=[1.0e8, 1.0e12], psym=-6, yrange=[-0.2,1.0]
  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 
  norm = 0.0 

  for i = 0, array_length-1 do begin 

     norm = norm + foo[i] 

     dndx[i] = foo[i] / dp3[i] 

  endfor 

  print, 'norm=', norm 

  m3 = m2 
  cbyfe3 = cbyfe2 
  dndx3 = dndx


; plot, cbyfe2, dndx, psym=-6

  mean = mean/ldla

  print, 'ldla=', ldla 
  print, 'mean=', mean

;--------------------------------------

  c_data2 = read_ascii('set49/halos_Zn.out', template=stars_template)
  fe_data2 = read_ascii('set49/halos_fe.out', template=stars_template)
  m_data2 = read_ascii('set49/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set49/halos_gas.out', template=stars_template)

  c2 = c_data2.field001[lcolumn:*,row]
  fe2 = fe_data2.field001[lcolumn:*,row]
  m2 = m_data2.field001[lcolumn:*,row]
  n2 = n_data2.field001[lcolumn:*,row]
  gas2 = gas_data2.field001[lcolumn:*,row]

  array_size = size(c2)
  array_length = array_size[1] 

  cbyfe2 = fltarr(array_length)
  febyh2 = fltarr(array_length)
  cross_section = fltarr(array_length)
  amplitude = fltarr(array_length)
  dndx = fltarr(array_length)
  foo = fltarr(array_length) 

;  sigma0 = 40.0                 ; kpc^2
;  m0 = 10.0^(-0.5)              ; 10^10 M_solar 
;  alpha = 0.2                   ; dimensionless 

  for i = 0, array_length-1 do begin 

     cnst = -0.05 

     if (i eq 0) then begin 
        dm = 0.01 * m2[i]
     endif else begin 
        dm = m2[i]-m2[i-1]
     endelse

     mH2 = gas2[i]*0.71
     
     if (fe2[i] eq 0.0) then begin 
        febyh2[i] = -10.0
     endif else begin 
        febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
     endelse
     
     if (c2[i] eq 0.0) then begin 
        cbyfe2[i] = -10.0
     endif else begin 
        cbyfe2[i] = alog10(abs(c2[i]/fe2[i])) + 3.07
     endelse

                                ; We need to find mass of a halo at
                                ; z=3 that has the same circular
                                ; velocity as this halo. 
     omega_mz = omega_nr*(1.0+z)^3/(omega_nr*(1.0+z)^3+omega_lambda+omega_k*(1.0+z)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     vc = 23.4 * (m2[i]*smallh/1.0e-2)^(1.0/3.0) $
          * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) $
          * ((1.0+z)/10.0)^0.5  ; km/s 

     omega_mz = omega_nr*(4.0)^3/(omega_nr*(4.0)^3+omega_lambda+omega_k*(4.0)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     m_z3 = vc/(((4.0)/10.0)^0.5 * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) * (smallh/1.0e-2)^(1.0/3.0) * 23.4)
     m_z3 = m_z3^3 

                                ; cross_section[i] = sigma0 * (m2[i]/m0)^2 * (1.0 + m2[i]/m0)^(alpha-2.0) ; kpc^2 
     cross_section[i] = sigma0 * (m_z3/m0)^2 * (1.0 + m_z3/m0)^(alpha-2.0) ; kpc^2 
     amplitude[i] = cross_section[i]*n2[i]*1.0e-6                          ; dimensionless (1.0e-6 is kpc^2/Mpc^2)

                                ; dndx is a misnomer here! What I mean
                                ; is (d^2n/dXdp . dp).  See discussion
                                ; around Eqn. 3 of Pontzen et al. (2008)
                                ; dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/(hubble_0*(1.0+z)^3) ; dimensionless 
     dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/hubble_0 ; dimensionless 
     foo[i] = dndx[i] 

  endfor

  deepee = fltarr(array_length) 

  mean = 0.0 
  ldla = 0.0 
  for i = 0, array_length-1 do begin 

     if (i eq 0) then begin 
        dp = abs(cbyfe2[i])
     endif else begin 
        dp = abs(cbyfe2[i] - cbyfe2[i-1])
     endelse

     deepee[i] = dp 

     if (dp EQ 0.0) then dp = 0.001

     dndx[i] = dndx[i]/dp

     ldla = ldla + (dndx[i] * dp) 

     mean = mean + (dndx[i] * cbyfe2[i] * dp)

  endfor 

  a = float(ylbound)
  b = float(yubound) 

  m2 = m2*1.0e10 
; plot, m2, cbyfe2, /xlog, xrange=[1.0e8, 1.0e12], psym=-6, yrange=[-0.2,1.0]
  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 
  norm = 0.0 

  for i = 0, array_length-1 do begin 

     norm = norm + foo[i] 

     dndx[i] = foo[i] / dp3[i] 

  endfor 

  print, 'norm=', norm 

  m1 = m2 
  cbyfe1 = cbyfe2 
  dndx1 = dndx


; plot, cbyfe2, dndx, psym=-6

  mean = mean/ldla

  print, 'ldla=', ldla 
  print, 'mean=', mean


;----------------------------------------

  c_data2 = read_ascii('set50/halos_Zn.out', template=stars_template)
  fe_data2 = read_ascii('set50/halos_fe.out', template=stars_template)
  m_data2 = read_ascii('set50/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set50/halos_gas.out', template=stars_template)

  c2 = c_data2.field001[lcolumn:*,row]
  fe2 = fe_data2.field001[lcolumn:*,row]
  m2 = m_data2.field001[lcolumn:*,row]
  n2 = n_data2.field001[lcolumn:*,row]
  gas2 = gas_data2.field001[lcolumn:*,row]

  array_size = size(c2)
  array_length = array_size[1] 

  cbyfe2 = fltarr(array_length)
  febyh2 = fltarr(array_length)
  cross_section = fltarr(array_length)
  amplitude = fltarr(array_length)
  dndx = fltarr(array_length)
  foo = fltarr(array_length) 

;  sigma0 = 40.0                 ; kpc^2
;  m0 = 10.0^(-0.5)              ; 10^10 M_solar 
;  alpha = 0.2                   ; dimensionless 

  for i = 0, array_length-1 do begin 

     cnst = -0.05 

     if (i eq 0) then begin 
        dm = 0.01 * m2[i]
     endif else begin 
        dm = m2[i]-m2[i-1]
     endelse

     mH2 = gas2[i]*0.71
     
     if (fe2[i] eq 0.0) then begin 
        febyh2[i] = -10.0
     endif else begin 
        febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
     endelse
     
     if (c2[i] eq 0.0) then begin 
        cbyfe2[i] = -10.0
     endif else begin 
        cbyfe2[i] = alog10(abs(c2[i]/fe2[i])) + 3.07
     endelse

                                ; We need to find mass of a halo at
                                ; z=3 that has the same circular
                                ; velocity as this halo. 
     omega_mz = omega_nr*(1.0+z)^3/(omega_nr*(1.0+z)^3+omega_lambda+omega_k*(1.0+z)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     vc = 23.4 * (m2[i]*smallh/1.0e-2)^(1.0/3.0) $
          * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) $
          * ((1.0+z)/10.0)^0.5  ; km/s 

     omega_mz = omega_nr*(4.0)^3/(omega_nr*(4.0)^3+omega_lambda+omega_k*(4.0)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     m_z3 = vc/(((4.0)/10.0)^0.5 * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) * (smallh/1.0e-2)^(1.0/3.0) * 23.4)
     m_z3 = m_z3^3 

                                ; cross_section[i] = sigma0 * (m2[i]/m0)^2 * (1.0 + m2[i]/m0)^(alpha-2.0) ; kpc^2 
     cross_section[i] = sigma0 * (m_z3/m0)^2 * (1.0 + m_z3/m0)^(alpha-2.0) ; kpc^2 
     amplitude[i] = cross_section[i]*n2[i]*1.0e-6                          ; dimensionless (1.0e-6 is kpc^2/Mpc^2)

                                ; dndx is a misnomer here! What I mean
                                ; is (d^2n/dXdp . dp).  See discussion
                                ; around Eqn. 3 of Pontzen et al. (2008)
                                ; dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/(hubble_0*(1.0+z)^3) ; dimensionless 
     dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/hubble_0 ; dimensionless 
     foo[i] = dndx[i] 

  endfor

  deepee = fltarr(array_length) 

  mean = 0.0 
  ldla = 0.0 
  for i = 0, array_length-1 do begin 

     if (i eq 0) then begin 
        dp = abs(cbyfe2[i])
     endif else begin 
        dp = abs(cbyfe2[i] - cbyfe2[i-1])
     endelse

     deepee[i] = dp 

     if (dp EQ 0.0) then dp = 0.001

     dndx[i] = dndx[i]/dp

     ldla = ldla + (dndx[i] * dp) 

     mean = mean + (dndx[i] * cbyfe2[i] * dp)

  endfor 

  a = float(ylbound)
  b = float(yubound) 

  m2 = m2*1.0e10 

  ms = rebin(m2, array_length/2)
  cs = rebin(cbyfe2, array_length/2)
  ms1 = rebin(m1, array_length/2)
  cs1 = rebin(cbyfe1, array_length/2)
  ms3 = rebin(m3, array_length/2)
  cs3 = rebin(cbyfe3, array_length/2)

;; plot, ms3, cs3, /xlog, xrange=[1.0e8, 1.0e12], psym=-6, yrange=[-1.0,0.0], xtitle='halo mass', ytitle='[Zn/Fe]', symsize=0.5
;; oplot, ms1, cs1, color=2, psym=-6, symsize=0.5
;; oplot, ms, cs, color=3, psym=-6, symsize=0.5

  ms_znbyfe = ms
  cs_znbyfe = cs
  ms1_znbyfe = ms1
  cs1_znbyfe = cs1
  ms3_znbyfe = ms3
  cs3_znbyfe = cs3

  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 

  norm = 0.0 
  for i = 0, array_length-1 do begin 
     norm = norm + foo[i] 
     dndx[i] = foo[i] / dp3[i] 
  endfor 
  print, 'norm=', norm 

  n = 50 
  t = fltarr(n)
  dx = (max(cbyfe2)-min(cbyfe2))/n 
  xmin = min(cbyfe2) 
  for i = 0, n-1 do begin 
     t[i] = xmin + dx*i 
  endfor 

  ns = rebin(dndx, array_length/2) 
  cs = rebin(cbyfe2, array_length/2)

;plot, cbyfe3, dndx3, psym=-6, xtitle='[Zn/Fe]', ytitle='d!E2!NN/dXd[Zn/Fe]', xrange=[-1.0,0.0], symsize=0.5


  ns1 = rebin(dndx1, array_length/2) 
  cs1 = rebin(cbyfe1, array_length/2)

;oplot, cs1, ns1, psym=-6, color=2, symsize=0.5

;oplot, cs, ns, psym=6, color=3, symsize=0.5

  ns_znbyfe = ns
  cs_znbyfe2 = cs
  ns1_znbyfe = ns1
  cs1_znbyfe2 = cs1
  ns3_znbyfe = dndx3
  cs3_znbyfe2 = cbyfe3

  mean = mean/ldla
  print, 'ldla=', ldla 
  print, 'mean=', mean


;; >>>>>3<<<<<

;------------------------------------------------------------------------

; IMPORTANT below:
; C = N
; Fe = O

  c_data2 = read_ascii('set95/halos_N.out', template=stars_template)
  fe_data2 = read_ascii('set95/halos_o.out', template=stars_template)
  m_data2 = read_ascii('set95/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set95/halos_gas.out', template=stars_template)

  c2 = c_data2.field001[lcolumn:*,row]
  fe2 = fe_data2.field001[lcolumn:*,row]
  m2 = m_data2.field001[lcolumn:*,row]
  n2 = n_data2.field001[lcolumn:*,row]
  gas2 = gas_data2.field001[lcolumn:*,row]

  array_size = size(c2)
  array_length = array_size[1] 

  cbyfe2 = fltarr(array_length)
  febyh2 = fltarr(array_length)
  cross_section = fltarr(array_length)
  amplitude = fltarr(array_length)
  dndx = fltarr(array_length)
  foo = fltarr(array_length) 

;  sigma0 = 40.0                 ; kpc^2
;  m0 = 10.0^(-0.5)              ; 10^10 M_solar 
;  alpha = 0.2                   ; dimensionless 

  for i = 0, array_length-1 do begin 

     cnst = -0.05 

     if (i eq 0) then begin 
        dm = 0.01 * m2[i]
     endif else begin 
        dm = m2[i]-m2[i-1]
     endelse

     mH2 = gas2[i]*0.71
     
     if (fe2[i] eq 0.0) then begin 
        febyh2[i] = -10.0
     endif else begin 
        febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
     endelse
     
     if (c2[i] eq 0.0) then begin 
        cbyfe2[i] = -10.0
     endif else begin 
        cbyfe2[i] = alog10(abs(c2[i]/fe2[i])) + 0.94
     endelse

                                ; We need to find mass of a halo at
                                ; z=3 that has the same circular
                                ; velocity as this halo. 
     omega_mz = omega_nr*(1.0+z)^3/(omega_nr*(1.0+z)^3+omega_lambda+omega_k*(1.0+z)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     vc = 23.4 * (m2[i]*smallh/1.0e-2)^(1.0/3.0) $
          * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) $
          * ((1.0+z)/10.0)^0.5  ; km/s 

     omega_mz = omega_nr*(4.0)^3/(omega_nr*(4.0)^3+omega_lambda+omega_k*(4.0)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     m_z3 = vc/(((4.0)/10.0)^0.5 * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) * (smallh/1.0e-2)^(1.0/3.0) * 23.4)
     m_z3 = m_z3^3 

                                ; cross_section[i] = sigma0 * (m2[i]/m0)^2 * (1.0 + m2[i]/m0)^(alpha-2.0) ; kpc^2 
     cross_section[i] = sigma0 * (m_z3/m0)^2 * (1.0 + m_z3/m0)^(alpha-2.0) ; kpc^2 
     amplitude[i] = cross_section[i]*n2[i]*1.0e-6                          ; dimensionless (1.0e-6 is kpc^2/Mpc^2)

                                ; dndx is a misnomer here! What I mean
                                ; is (d^2n/dXdp . dp).  See discussion
                                ; around Eqn. 3 of Pontzen et al. (2008)
                                ; dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/(hubble_0*(1.0+z)^3) ; dimensionless 
     dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/hubble_0 ; dimensionless 
     foo[i] = dndx[i] 

  endfor

  deepee = fltarr(array_length) 

  mean = 0.0 
  ldla = 0.0 
  for i = 0, array_length-1 do begin 

     if (i eq 0) then begin 
        dp = abs(cbyfe2[i])
     endif else begin 
        dp = abs(cbyfe2[i] - cbyfe2[i-1])
     endelse

     deepee[i] = dp 

     if (dp EQ 0.0) then dp = 0.001

     dndx[i] = dndx[i]/dp

     ldla = ldla + (dndx[i] * dp) 

     mean = mean + (dndx[i] * cbyfe2[i] * dp)

  endfor 

  a = float(ylbound)
  b = float(yubound) 

  m2 = m2*1.0e10 
; plot, m2, cbyfe2, /xlog, xrange=[1.0e8, 1.0e12], psym=-6, yrange=[-0.2,1.0]
  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 
  norm = 0.0 

  for i = 0, array_length-1 do begin 

     norm = norm + foo[i] 

     dndx[i] = foo[i] / dp3[i] 

  endfor 

  print, 'norm=', norm 

  m3 = m2 
  cbyfe3 = cbyfe2 
  dndx3 = dndx


; plot, cbyfe2, dndx, psym=-6

  mean = mean/ldla

  print, 'ldla=', ldla 
  print, 'mean=', mean

;--------------------------------------

  c_data2 = read_ascii('set49/halos_N.out', template=stars_template)
  fe_data2 = read_ascii('set49/halos_o.out', template=stars_template)
  m_data2 = read_ascii('set49/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set49/halos_gas.out', template=stars_template)

  c2 = c_data2.field001[lcolumn:*,row]
  fe2 = fe_data2.field001[lcolumn:*,row]
  m2 = m_data2.field001[lcolumn:*,row]
  n2 = n_data2.field001[lcolumn:*,row]
  gas2 = gas_data2.field001[lcolumn:*,row]

  array_size = size(c2)
  array_length = array_size[1] 

  cbyfe2 = fltarr(array_length)
  febyh2 = fltarr(array_length)
  cross_section = fltarr(array_length)
  amplitude = fltarr(array_length)
  dndx = fltarr(array_length)
  foo = fltarr(array_length) 

;  sigma0 = 40.0                 ; kpc^2
;  m0 = 10.0^(-0.5)              ; 10^10 M_solar 
;  alpha = 0.2                   ; dimensionless 

  for i = 0, array_length-1 do begin 

     cnst = -0.05 

     if (i eq 0) then begin 
        dm = 0.01 * m2[i]
     endif else begin 
        dm = m2[i]-m2[i-1]
     endelse

     mH2 = gas2[i]*0.71
     
     if (fe2[i] eq 0.0) then begin 
        febyh2[i] = -10.0
     endif else begin 
        febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
     endelse
     
     if (c2[i] eq 0.0) then begin 
        cbyfe2[i] = -10.0
     endif else begin 
        cbyfe2[i] = alog10(abs(c2[i]/fe2[i])) + 0.94
     endelse

                                ; We need to find mass of a halo at
                                ; z=3 that has the same circular
                                ; velocity as this halo. 
     omega_mz = omega_nr*(1.0+z)^3/(omega_nr*(1.0+z)^3+omega_lambda+omega_k*(1.0+z)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     vc = 23.4 * (m2[i]*smallh/1.0e-2)^(1.0/3.0) $
          * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) $
          * ((1.0+z)/10.0)^0.5  ; km/s 

     omega_mz = omega_nr*(4.0)^3/(omega_nr*(4.0)^3+omega_lambda+omega_k*(4.0)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     m_z3 = vc/(((4.0)/10.0)^0.5 * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) * (smallh/1.0e-2)^(1.0/3.0) * 23.4)
     m_z3 = m_z3^3 

                                ; cross_section[i] = sigma0 * (m2[i]/m0)^2 * (1.0 + m2[i]/m0)^(alpha-2.0) ; kpc^2 
     cross_section[i] = sigma0 * (m_z3/m0)^2 * (1.0 + m_z3/m0)^(alpha-2.0) ; kpc^2 
     amplitude[i] = cross_section[i]*n2[i]*1.0e-6                          ; dimensionless (1.0e-6 is kpc^2/Mpc^2)

                                ; dndx is a misnomer here! What I mean
                                ; is (d^2n/dXdp . dp).  See discussion
                                ; around Eqn. 3 of Pontzen et al. (2008)
                                ; dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/(hubble_0*(1.0+z)^3) ; dimensionless 
     dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/hubble_0 ; dimensionless 
     foo[i] = dndx[i] 

  endfor

  deepee = fltarr(array_length) 

  mean = 0.0 
  ldla = 0.0 
  for i = 0, array_length-1 do begin 

     if (i eq 0) then begin 
        dp = abs(cbyfe2[i])
     endif else begin 
        dp = abs(cbyfe2[i] - cbyfe2[i-1])
     endelse

     deepee[i] = dp 

     if (dp EQ 0.0) then dp = 0.001

     dndx[i] = dndx[i]/dp

     ldla = ldla + (dndx[i] * dp) 

     mean = mean + (dndx[i] * cbyfe2[i] * dp)

  endfor 

  a = float(ylbound)
  b = float(yubound) 

  m2 = m2*1.0e10 
; plot, m2, cbyfe2, /xlog, xrange=[1.0e8, 1.0e12], psym=-6, yrange=[-0.2,1.0]
  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 
  norm = 0.0 

  for i = 0, array_length-1 do begin 

     norm = norm + foo[i] 

     dndx[i] = foo[i] / dp3[i] 

  endfor 

  print, 'norm=', norm 

  m1 = m2 
  cbyfe1 = cbyfe2 
  dndx1 = dndx


; plot, cbyfe2, dndx, psym=-6

  mean = mean/ldla

  print, 'ldla=', ldla 
  print, 'mean=', mean


;----------------------------------------

  c_data2 = read_ascii('set50/halos_N.out', template=stars_template)
  fe_data2 = read_ascii('set50/halos_o.out', template=stars_template)
  m_data2 = read_ascii('set50/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set50/halos_gas.out', template=stars_template)

  c2 = c_data2.field001[lcolumn:*,row]
  fe2 = fe_data2.field001[lcolumn:*,row]
  m2 = m_data2.field001[lcolumn:*,row]
  n2 = n_data2.field001[lcolumn:*,row]
  gas2 = gas_data2.field001[lcolumn:*,row]

  array_size = size(c2)
  array_length = array_size[1] 

  cbyfe2 = fltarr(array_length)
  febyh2 = fltarr(array_length)
  cross_section = fltarr(array_length)
  amplitude = fltarr(array_length)
  dndx = fltarr(array_length)
  foo = fltarr(array_length) 

;  sigma0 = 40.0                 ; kpc^2
;  m0 = 10.0^(-0.5)              ; 10^10 M_solar 
;  alpha = 0.2                   ; dimensionless 

  for i = 0, array_length-1 do begin 

     cnst = -0.05 

     if (i eq 0) then begin 
        dm = 0.01 * m2[i]
     endif else begin 
        dm = m2[i]-m2[i-1]
     endelse

     mH2 = gas2[i]*0.71
     
     if (fe2[i] eq 0.0) then begin 
        febyh2[i] = -10.0
     endif else begin 
        febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
     endelse
     
     if (c2[i] eq 0.0) then begin 
        cbyfe2[i] = -10.0
     endif else begin 
        cbyfe2[i] = alog10(abs(c2[i]/fe2[i])) + 0.94
     endelse

                                ; We need to find mass of a halo at
                                ; z=3 that has the same circular
                                ; velocity as this halo. 
     omega_mz = omega_nr*(1.0+z)^3/(omega_nr*(1.0+z)^3+omega_lambda+omega_k*(1.0+z)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     vc = 23.4 * (m2[i]*smallh/1.0e-2)^(1.0/3.0) $
          * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) $
          * ((1.0+z)/10.0)^0.5  ; km/s 

     omega_mz = omega_nr*(4.0)^3/(omega_nr*(4.0)^3+omega_lambda+omega_k*(4.0)^2) 
     d = omega_mz-1.0 
     delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
     m_z3 = vc/(((4.0)/10.0)^0.5 * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) * (smallh/1.0e-2)^(1.0/3.0) * 23.4)
     m_z3 = m_z3^3 

                                ; cross_section[i] = sigma0 * (m2[i]/m0)^2 * (1.0 + m2[i]/m0)^(alpha-2.0) ; kpc^2 
     cross_section[i] = sigma0 * (m_z3/m0)^2 * (1.0 + m_z3/m0)^(alpha-2.0) ; kpc^2 
     amplitude[i] = cross_section[i]*n2[i]*1.0e-6                          ; dimensionless (1.0e-6 is kpc^2/Mpc^2)

                                ; dndx is a misnomer here! What I mean
                                ; is (d^2n/dXdp . dp).  See discussion
                                ; around Eqn. 3 of Pontzen et al. (2008)
                                ; dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/(hubble_0*(1.0+z)^3) ; dimensionless 
     dndx[i] = amplitude[i]*speed_of_light*yrbys*cmbympc/hubble_0 ; dimensionless 
     foo[i] = dndx[i] 

  endfor

  deepee = fltarr(array_length) 

  mean = 0.0 
  ldla = 0.0 
  for i = 0, array_length-1 do begin 

     if (i eq 0) then begin 
        dp = abs(cbyfe2[i])
     endif else begin 
        dp = abs(cbyfe2[i] - cbyfe2[i-1])
     endelse

     deepee[i] = dp 

     if (dp EQ 0.0) then dp = 0.001

     dndx[i] = dndx[i]/dp

     ldla = ldla + (dndx[i] * dp) 

     mean = mean + (dndx[i] * cbyfe2[i] * dp)

  endfor 

  a = float(ylbound)
  b = float(yubound) 

  m2 = m2*1.0e10 

  ms = rebin(m2, array_length/2)
  cs = rebin(cbyfe2, array_length/2)
  ms1 = rebin(m1, array_length/2)
  cs1 = rebin(cbyfe1, array_length/2)
  ms3 = rebin(m3, array_length/2)
  cs3 = rebin(cbyfe3, array_length/2)

;; plot, ms3, cs3, /xlog, xrange=[1.0e8, 1.0e12], psym=-6, yrange=[-3.0,-2.6], xtitle='halo mass', ytitle='[N/O]', symsize=0.5
;; oplot, ms1, cs1, color=2, psym=-6, symsize=0.5
;; oplot, ms, cs, color=3, psym=-6, symsize=0.5

  ms_nbyo = ms
  cs_nbyo = cs
  ms1_nbyo = ms1
  cs1_nbyo = cs1
  ms3_nbyo = ms3
  cs3_nbyo = cs3

  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 

  norm = 0.0 
  for i = 0, array_length-1 do begin 
     norm = norm + foo[i] 
     dndx[i] = foo[i] / dp3[i] 
  endfor 
  print, 'norm=', norm 

  n = 50 
  t = fltarr(n)
  dx = (max(cbyfe2)-min(cbyfe2))/n 
  xmin = min(cbyfe2) 
  for i = 0, n-1 do begin 
     t[i] = xmin + dx*i 
  endfor 


  ns = rebin(dndx, array_length/2) 
  cs = rebin(cbyfe2, array_length/2)

; plot, cbyfe3, dndx3, psym=-6, xtitle='[N/O]', ytitle='d!E2!NN/dXd[N/O]', xrange=[-3.0,-2.6], symsize=0.5



  ns1 = rebin(dndx1, array_length/2) 
  cs1 = rebin(cbyfe1, array_length/2)

; oplot, cs1, ns1, psym=-6, color=2, symsize=0.5



; oplot, cs, ns, psym=6, color=3, symsize=0.5
  mean = mean/ldla
  print, 'ldla=', ldla 
  print, 'mean=', mean

  ns_nbyo = ns
  cs_nbyo2 = cs
  ns1_nbyo = ns1
  cs1_nbyo2 = cs1
  ns3_nbyo = dndx3
  cs3_nbyo2 = cbyfe3



; XYOuts, 0.5, 0.95, ALIGNMENT=0.5, CHARSIZE=2.0, /NORMAL, redshift_label

;; >>>>>4<<<<<

;------------------------------------------------------------------------------

; Plot everything

  plot, ms3_cbyfe, cs3_cbyfe, /xlog, xrange=[1.0e8, 1.0e11], yrange=[-0.2,1.0], xtitle='!6halo mass', ytitle='[C/Fe]', symsize=0.5, thick=3
  oplot, ms1_cbyfe, cs1_cbyfe, color=2, symsize=0.5, thick=3
  oplot, ms_cbyfe, cs_cbyfe, color=3, symsize=0.5, thick=3
  xyouts, 1.0e10, 0.8, redshift_label, charsize=0.7

  plot, ms3_obysi, cs3_obysi, /xlog, xrange=[1.0e8, 1.0e11], yrange=[-0.5,1.0], xtitle='halo mass', ytitle='[O/Si]', symsize=0.5, thick=3
  oplot, ms1_obysi, cs1_obysi, color=2, symsize=0.5, thick=3
  oplot, ms_obysi, cs_obysi, color=3, symsize=0.5, thick=3
  xyouts, 1.0e10, 0.75, redshift_label, charsize=0.7

  plot, ms3_znbyfe, cs3_znbyfe, /xlog, xrange=[1.0e8, 1.0e11], yrange=[-1.0,0.0], xtitle='halo mass', ytitle='[Zn/Fe]', symsize=0.5, thick=3
  oplot, ms1_znbyfe, cs1_znbyfe, color=2, symsize=0.5, thick=3
  oplot, ms_znbyfe, cs_znbyfe, color=3, symsize=0.5, thick=3
  xyouts, 1.0e10, -0.8, redshift_label, charsize=0.7

  plot, ms3_nbyo, cs3_nbyo, /xlog, xrange=[1.0e8, 1.0e11], yrange=[-3.0,-2.6], xtitle='halo mass', ytitle='[N/O]', symsize=0.5, thick=3
  oplot, ms1_nbyo, cs1_nbyo, color=2, symsize=0.5, thick=3
  oplot, ms_nbyo, cs_nbyo, color=3, symsize=0.5, thick=3
  xyouts, 1.0e10, -2.95, redshift_label, charsize=0.7

  plot, cs3_cbyfe2, ns3_cbyfe, xtitle='[C/Fe]', ytitle='d!E2!NN/dXd[C/Fe]', xrange=[-0.2,1.0], symsize=0.5, thick=3
  oplot, cs1_cbyfe2, ns1_cbyfe, color=2, symsize=0.5, thick=3 
  oplot, cs_cbyfe2, ns_cbyfe, color=3, symsize=0.5, thick=3
; xyouts, 0.8, 0.006, redshift_label, charsize=0.7
  xyouts, 0.8, 0.5, redshift_label, charsize=0.7

  plot, cs3_obysi2, ns3_obysi, xtitle='[O/Si]', ytitle='d!E2!NN/dXd[O/Si]', xrange=[-0.5,1.0], symsize=0.5, thick=3
  oplot, cs1_obysi2, ns1_obysi, color=2, symsize=0.5, thick=3
  oplot, cs_obysi2, ns_obysi, color=3, symsize=0.5, thick=3

  array_size = size(cs_obysi2)
  array_length = array_size[1] 
;  print, array_length 

;  for i = 0, array_length-1 do begin 
;     print, cs_obysi2[i], ns_obysi[i]
;  endfor 

;xyouts, 0.75, 0.006, redshift_label, charsize=0.7
  xyouts, 0.75, 0.5, redshift_label, charsize=0.7

  plot, cs3_znbyfe2, ns3_znbyfe, xtitle='[Zn/Fe]', ytitle='d!E2!NN/dXd[Zn/Fe]', xrange=[-1.0,0.0], symsize=0.5, thick=3
  oplot, cs1_znbyfe2, ns1_znbyfe, color=2, symsize=0.5, thick=3
  oplot, cs_znbyfe2, ns_znbyfe, color=3, symsize=0.5, thick=3
;xyouts, -0.8, 0.006, redshift_label, charsize=0.7
  xyouts, -0.8, 0.5, redshift_label, charsize=0.7

  plot, cs3_nbyo2, ns3_nbyo, xtitle='[N/O]', ytitle='d!E2!NN/dXd[N/O]', xrange=[-3.0,-2.6], symsize=0.5, thick=3
  oplot, cs1_nbyo2, ns1_nbyo, color=2, symsize=0.5, thick=3
  oplot, cs_nbyo2, ns_nbyo, color=3, symsize=0.5, thick=3
;xyouts, -2.93, 0.006, redshift_label, charsize=0.7
  xyouts, -2.93, 0.5, redshift_label, charsize=0.7

;; device, /close_file
;; set_plot, 'X'

END
