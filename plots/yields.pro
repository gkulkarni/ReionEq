
set_plot, 'ps'
device, filename='yields.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

;; window, xsize=1000, ysize=1000
;; Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 128, 128, 4
TvLCT, 0, 255, 0, 5
TvLCT, 255, 128, 0, 6
TvLCT, 128, 0, 128, 7
TvLCT, 255, 0, 255, 8
TvLCT, 255, 255, 0, 9
!P.charsize = 1.5
!P.thick=2

readcol, '../data/yields.t.wwhw.Z0', m, mrem, ekin, vinf, d, he, $
         c, n, o, mg, si, s, ca, fe, zn, z 

plot, m, c, psym=-6, /xlog, /ylog, /nodata, yrange=[1.0e-30,1.e2], $
      xtitle='!6initial stellar mass (M!D!9n!X!N)', ytitle='!6metal ejecta mass (M!D!9n!X!N)'
oplot, m, c, psym=-6, color=2
oplot, m, n, psym=-6, color=3
oplot, m, o, psym=-6, color=4
oplot, m, si, psym=-6, color=5
oplot, m, fe, psym=-6, color=6
oplot, m, zn, psym=-6, color=7
legend, ['C', 'N', 'O', 'Si', 'Fe', 'Zn'], linestyle=[0,0,0,0,0,0], $
        psym=[-6,-6,-6,-6,-6,-6], color=[2,3,4,5,6,7], /bottom, /right

;; c = c/m 
;; n = n/m 
;; o = o/m 
;; si = si/m 
;; fe = fe/m 
;; zn = zn/m 

;; plot, m, c, psym=-6, /xlog, /ylog, /nodata, yrange=[1.0e-30,1.e2]
;; oplot, m, c, psym=-6, color=2
;; oplot, m, n, psym=-6, color=3
;; oplot, m, o, psym=-6, color=4
;; oplot, m, si, psym=-6, color=5
;; oplot, m, fe, psym=-6, color=6
;; oplot, m, zn, psym=-6, color=7
;; legend, ['C', 'N', 'O', 'Si', 'Fe', 'Zn'], linestyle=[0,0,0,0,0,0], psym=[-6,-6,-6,-6,-6,-6], color=[2,3,4,5,6,7], charsize=1.2, /bottom, /right

device, /close_file
set_plot, 'X'


