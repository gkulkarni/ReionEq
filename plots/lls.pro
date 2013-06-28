window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4 
!P.charsize = 2

readcol, 'set68/reion.out', z, q, tau, gammapi, temph, $
         tempc, avtemp, x_ii, dnlldz, $
         lmfp, r, igmdcrit, nphdot, temphva, fv, /silent 
plot, z, dnlldz, xrange=[2,16], xstyle=1, /ylog

readcol, 'WMAP5atomic_IGM.dat', z, q, gpi1, gpi2, gpi3, t1, $
         t2, cp, t3, t4, dn, mfp, tau, xglob, ngamma 
oplot, z, dn, color=4

;; readcol, '10^10_IGM.dat', z, q, gpi1, gpi2, gpi3, t1, $
;;          t2, cp, t3, t4, dn, mfp, tau, xhi, ngamma, t5
;; oplot, z, dn, color=4




