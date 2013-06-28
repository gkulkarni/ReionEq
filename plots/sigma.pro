
; File: sigma.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Computes the HI absorption cross-section of a given halo.  Is used
; in zism_cs.pro.

FUNCTION sigma, m, z  

  ; m = halo mass 
  ; z = halo redshift 
  ; output sigma = cross section in kpc^2 

  smallh = 0.71 
  yrbys = 3.154e7
  cmbympc = 3.241e-25 
  hubble_0 = 1.023e-10*0.71     ; yr^-1 
  omega_nr = 0.2646
  omega_lambda = 0.734
  omega_k = 1.0-(omega_lambda+omega_nr)  
  pi = 3.14159265358979
  speed_of_light = 2.998e10     ; cm/s
  sigma0 = 40.0                 ; kpc^2
  m0 = 10.0^(-0.5)              ; 10^10 M_solar 
  alpha = 0.5                   ; dimensionless 


                                ; We need to find mass of a halo at
                                ; z=3 that has the same circular
                                ; velocity as this halo. 
  omega_mz = omega_nr*(1.0+z)^3/(omega_nr*(1.0+z)^3+omega_lambda+omega_k*(1.0+z)^2) 
  d = omega_mz-1.0 
  delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
  vc = 23.4 * (m*smallh/1.0e-2)^(1.0/3.0) $
       * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) $
       * ((1.0+z)/10.0)^0.5     ; km/s 
  omega_mz = omega_nr*(4.0)^3/(omega_nr*(4.0)^3+omega_lambda+omega_k*(4.0)^2) 
  d = omega_mz-1.0 
  delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d 
  m_z3 = vc/(((4.0)/10.0)^0.5 * ((omega_nr/omega_mz)*(delta_c/(18.0*pi*pi)))^(1.0/6.0) * (smallh/1.0e-2)^(1.0/3.0) * 23.4)
  m_z3 = m_z3^3 

  sigma = sigma0 * (m_z3/m0)^2 * (1.0 + m_z3/m0)^(alpha-2.0) ; kpc^2 
  return, sigma 

END 
