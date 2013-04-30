; This code snippet was used in lf.pro to understand AB magnitudes.

  dl = 106251.8d0 ; Luminosity distance to z=10 (Mpc) 
  cmbympc = 3.24077928965d-25
  dl_cm = dl / cmbympc ; Luminosity distance to z=10 (cm) 
  pi = 3.14159265358979d0 
  dl_pc = dl * 1.0e6 ; Luminosity distance to z=10 (pc) 

  for i = 0, total_indices-1 do begin 
     f = 10.0d0^((51.6d0-l[i])/2.5d0) ; erg/s/Hz
     l_alt = -2.5d0*alog10(f/(4.0d0*pi*dl_cm*dl_cm))-48.60d0
     l_altabs = l_alt - 5.0d0*alog10(dl_pc/10.0d0)
     chk = 2.5d0*alog10(4.0d0*pi) + 5.0d0*alog10(1.0/(cmbympc*1.0d6)) - 48.60d0 + 5.0*alog10(10.0d0) 
     print, i, l[i], f, l_alt, l_altabs, chk 
  endfor

