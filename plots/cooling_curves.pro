
PRO cooling_curves 

  ; Plot various cooling curves as a check.

  window, xsize=1000, ysize=1000
  TvLCT, 255, 0, 0, 2 
  Device, decomposed=0
  !P.charsize = 2.0

  jbyerg = 1.0d7
  mbycm_cubed = 1.0e6

  readcol, '../data/SD93-00.cie', t, lambda, format='(D,D)', /silent 
  lambda = lambda*jbyerg*mbycm_cubed 
  plot, t, lambda, /xlog, /ylog, xtitle='T (K)', ytitle='!7K!X (erg cm!E3!N s!E-1!N)', $
        yrange=[1.0d-25, 1.0d-19]

  readcol, '../data/SD93+05.cie', t, lambda, format='(D,D)', /silent 
  lambda = lambda*jbyerg*mbycm_cubed 
  oplot, t, lambda

  readcol, '../data/SD93-10.cie', t, lambda, format='(D,D)', /silent 
  lambda = lambda*jbyerg*mbycm_cubed 
  oplot, t, lambda

  readcol, '../data/SD93-15.cie', t, lambda, format='(D,D)', /silent 
  lambda = lambda*jbyerg*mbycm_cubed 
  oplot, t, lambda

  readcol, '../data/SD93-20.cie', t, lambda, format='(D,D)', /silent 
  lambda = lambda*jbyerg*mbycm_cubed 
  oplot, t, lambda

  readcol, '../data/SD93-30.cie', t, lambda, format='(D,D)', /silent 
  lambda = lambda*jbyerg*mbycm_cubed 
  oplot, t, lambda

  readcol, '../data/SD93_zero.cie', t, lambda, format='(D,D)', /silent 
  lambda = lambda*jbyerg*mbycm_cubed 
  oplot, t, lambda

  legend, ['(each curve corresponds to an [Fe/H];', 'all curves from Sutherland and Dopita 93)'], $
          /bottom, /right, box=0

  legend, ['collisional ionization equilibrium (CIE)', 'non-equilibrium (NEQ)'], $
          linestyle=[0,0], color=[-1,2], /right, box=0

  readcol, '../data/SD93_nie_zero.dat', t, lambda, format='(D,D)', /silent 
  lambda = lambda*jbyerg*mbycm_cubed 
  oplot, t, lambda, color=2
  
  readcol, '../data/SD93_nie_metals.dat', t, lambda, format='(D,D)', /silent 
  lambda = lambda*jbyerg*mbycm_cubed 
  oplot, t, lambda, color=2

END 

