
readcol, 'set68/reion.out', z, q, tau, gammapi, temph, $
         tempc, (q*temph+(1.0_prec-q)*tempc), x_ii, dnlldz, $
         lmfp, r, igmdcrit, nphdot, temphva, fv, /silent 
plot, z, dnlldz 
