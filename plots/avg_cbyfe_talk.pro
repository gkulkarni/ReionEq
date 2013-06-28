set_plot, 'ps'
device, filename='avg_cbyfe_talk.ps', xsize=7.0, ysize=3.0, $
        /inches, color=1, /HELVETICA
!P.font=0

;; window, xsize=500, ysize=300
;; Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
!P.charsize = 1.5
!P.thick = 3.0
!P.charthick = 1

rs = fltarr(10)
cbyfe_avg = fltarr(10)
plot, rs, cbyfe_avg, yrange=[-1.0,1.0], ytitle='[C/Fe]', $
      xrange=[2,6.5], xstyle=1, /nodata, xthick=3, ythick=3 

readcol, 'avg_data/cs2_51_cbyfe.dat', rs, cbyfe_avg, /silent 
oplot, rs, cbyfe_avg, thick=6
readcol, 'avg_data/cs2_49_cbyfe.dat', rs, cbyfe_avg, /silent 
oplot, rs, cbyfe_avg, color=2, thick=6
readcol, 'avg_data/cs2_50_cbyfe.dat', rs, cbyfe_avg, /silent 
oplot, rs, cbyfe_avg, color=3, thick=6

plotsym, 0, 0.7, /FILL
readcol, '../data/abr_cooke_cbyfe.dat', z, c, cl, cu 
oplot, z, c, psym=8 
dy1 = cu - c 
dy2 = c - cl 
oploterror, z, c, dy1, psym=8, /hibar
oploterror, z, c, dy2, psym=8, /lobar
readcol, '../data/abr_becker_cbyfe.dat', z, c, cerr 
zt = fltarr(4) 
ct = fltarr(4) 
zt2 = zt 
ct2 = ct 
cerrt = fltarr(4) 
j = 0 
k = 0 
for i = 0, 7 do begin 
   if cerr[i] ne 0.0 then begin 
      zt[j] = z[i] 
      ct[j] = c[i] 
      cerrt[j] = cerr[i] 
      j = j+1 
   endif else begin 
      zt2[k] = z[i] 
      ct2[k] = c[i] 
      k = k+1 
   endelse
endfor
oploterror, zt, ct, cerrt, psym=8, thick=3
plotsym, 2, 2
oplot, zt2, ct2, psym=8, thick=3

device, /close_file
set_plot, 'X'

END

