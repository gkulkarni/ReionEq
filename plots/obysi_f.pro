; IMPORTANT below:
; C = O 
; Fe = Si 

  c_data2 = read_ascii('set125/halos_o.out', template=stars_template)
  fe_data2 = read_ascii('set125/halos_Si.out', template=stars_template)
  m_data2 = read_ascii('set125/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set125/halos_gas.out', template=stars_template)

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
  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 
  norm = 0.0 

  for i = 0, array_length-1 do begin 

     norm = norm + foo[i] 

     dndx[i] = foo[i] / dp3[i] 

  endfor 

  m3 = m2 
  cbyfe3 = cbyfe2 
  dndx3 = dndx

;--------------------------------------

  c_data2 = read_ascii('set126/halos_o.out', template=stars_template)
  fe_data2 = read_ascii('set126/halos_Si.out', template=stars_template)
  m_data2 = read_ascii('set126/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set126/halos_gas.out', template=stars_template)

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
  dp2 = deepee
  dp2 = smooth(deepee, 80, /edge_truncate)
  dp3 = ts_smooth(deepee, 80, /forward) 
  norm = 0.0 

  for i = 0, array_length-1 do begin 

     norm = norm + foo[i] 

     dndx[i] = foo[i] / dp3[i] 

  endfor 

  m1 = m2 
  cbyfe1 = cbyfe2 
  dndx1 = dndx

;----------------------------------------

  ms1 = m1 ; rebin(m1, array_length/2)
  cs1 = cbyfe1 ; rebin(cbyfe1, array_length/2)
  ms3 = m3 ; rebin(m3, array_length/2)
  cs3 = cbyfe3 ; rebin(cbyfe3, array_length/2)

  ms1_obysi = ms1
  cs1_obysi = cs1
  ms3_obysi = ms3
  cs3_obysi = cs3

  ns1 = dndx1 ; rebin(dndx1, array_length/2) 
  cs1 = cbyfe1 ; rebin(cbyfe1, array_length/2)

  ns1_obysi = ns1
  cs1_obysi2 = cs1
  ns3_obysi = dndx3
  cs3_obysi2 = cbyfe3

