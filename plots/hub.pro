FUNCTION hub, z 
  
  smallh = 0.71 
  omega_nr = 0.2646
  omega_lambda = 0.734
  hubble_0 = 1.023e-10 * smallh
  hubp = hubble_0*sqrt(omega_nr*(1.0+z)^3+omega_lambda) ; yr^-1 
  return, hubp

END
