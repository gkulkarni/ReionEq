row = 87 
  z = 6.0

  ylbound = 1.0e-6
  yubound = 1.0e0

  window, xsize=1000, ysize=1000
  !P.multi=0
  Device, decomposed=0
  TvLCT, 255, 0, 0, 2 
  TvLCT, 0, 127, 255, 3
  !P.charsize = 1.5
  !P.thick = 1.0
  !P.charthick = 1

  redshift_label = 'z='+strtrim(uint(z),2)

  smallh = 0.71 
  lcolumn = 89
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

  c_data2 = read_ascii('set164/halos_o.out', template=stars_template)
  fe_data2 = read_ascii('set164/halos_Si.out', template=stars_template)
  m_data2 = read_ascii('set164/halos.out', template=stars_template) 
  n_data2 = read_ascii('set34/nofmh.out', template=stars_template) 
  gas_data2 = read_ascii('set164/halos_gas.out', template=stars_template)

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

;   plot, cbyfe2, dndx, psym=-6

  m3 = m2 
  cbyfe3 = cbyfe2 
  dndx3 = dndx

  ms3 = rebin(m3, array_length/2)
  cs3 = rebin(cbyfe3, array_length/2)
  ms3_obysi = ms3
  cs3_obysi = cs3
  ns3 = rebin(dndx3, array_length/2) 
  ; ns3_obysi = dndx3
  ns3_obysi = ns3
  ; cs3_obysi2 = cbyfe3
  cs3_obysi2 = cs3

  plot, ms3_obysi, cs3_obysi, /xlog, yrange=[-0.5,1.0]

   ; plot, cs3_obysi2, ns3_obysi, xtitle='[O/Si]', ytitle='d!E2!NN/dXd[O/Si]', xrange=[-0.5,0.5]

end
