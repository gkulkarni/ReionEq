; File: vir.pro 
;  Cre: 2013-03-05
;  Mod: $Date$ ($Revision$) 
; 
; Just a small IDL code to help understand the redshift evolution of
; v_circ, r_vir, and t_dyn. 


window, xsize=1000, ysize=1000
Device, decomposed=0
!P.charsize=2

omega_nr = 0.2646 
omega_lambda = 0.734
omega_k = 1.0 - (omega_lambda + omega_nr) 
pi = 3.14159265358979

zmin = 0.1
zmax = 100.0
incr = 1.02 

n = uint(alog10(zmax/zmin)/alog10(incr))
z = fltarr(n+1) 
vcirc = fltarr(n+1) 
rvir = fltarr(n+1)
tdyn = fltarr(n+1) 
tdyn2 = fltarr(n+1)

z[0] = zmin 

for i = 0, n do begin 

   omega_mz = omega_nr*(1.0+z[i])^3 / (omega_nr*(1.0+z[i])^3 + omega_lambda + omega_k*(1.0+z[i])^2) ; dimensionless 
   d = omega_mz - 1.0 ; dimensionless 
   delta_c = 18.0*pi*pi + 82.0*d + 39.0*d*d ; dimensionless

   vcirc[i] = (omega_nr*delta_c/(omega_mz*18.0*pi*pi))^(1.0/6.0) * (1.0+z[i])^0.5
   rvir[i] = (omega_nr*delta_c/(omega_mz*18.0*pi*pi))^(-1.0/3.0) * (1.0+z[i])^(-1.0)

   tdyn[i] = rvir[i]/vcirc[i]
   ; tdyn2[i] = rvir[i]/(vcirc[i]^2)

   tdyn2[i] = 1.0 / (1.0+2.0e4*exp(-z[i]^2/4.0))

   if i lt n then z[i+1] = z[i]*incr 

endfor

;; plot, z, vcirc, /xlog, /ylog
;; oplot, z, rvir, linestyle=5

plot, z, tdyn, /xlog, /ylog 
oplot, z, tdyn2, linestyle=5

END
