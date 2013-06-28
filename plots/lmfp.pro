window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4 
!P.charsize = 2

; Plot our result. 
readcol, 'set68/reion.out', z, q, tau, gammapi, temph, $
         tempc, avtemp, x_ii, dnlldz, $
         lmfp, r, igmdcrit, nphdot, temphva, fv, /silent 
plot, z, lmfp, /xlog, /ylog, xrange=[2,7], xstyle=1, $
      yrange=[5.0,300.0], xtitle='z', ytitle='mean free path (proper Mpc)' , ystyle=1

; Plot Tirth's result. 
;; readcol, 'WMAP5atomic_IGM.dat', z, q, gpi1, gpi2, gpi3, t1, $
;;          t2, cp, t3, t4, dn, mfp, tau, xglob, ngamma 
;; oplot, z, mfp, color=4, linestyle=1

;; readcol, '10^10_IGM.dat', z, q, gpi1, gpi2, gpi3, t1, $
;;          t2, cp, t3, t4, dn, mfp, tau, xhi, ngamma, t5
;; oplot, z, mfp, color=4

; Plot CAFG 2008 result. 
lmfp_fg = lmfp 
n = size(lmfp, /n_elements)
for i = 0, n-1 do begin 
   lmfp_fg[i] = 85.0*((1.0+z[i])/4.0)^(-4.0)
endfor
oplot, z, lmfp_fg, color=2, thick=3 

; Plot Madau 1999 result. 
lmfp_madau = lmfp 
for i = 0, n-1 do begin 
   lmfp_fg[i] = 33.0*((1.0+z[i])/4.0)^(-4.5)
endfor
oplot, z, lmfp_fg, color=3, thick=3

; Add legend. 
legend, ['my result', 'CAFG 08', 'Madau 99'], linestyle=[0,0,0], color=[-1,2,3], $
        thick=[1,3,3], /right

END
