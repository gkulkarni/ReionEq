
window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 127, 0, 153, 4
TvLCT, 255, 255, 0, 5
TvLCT, 1, 3, 64, 6
!P.charsize = 2.0
!P.thick = 1

readcol, '../dndm_PopIII_salpeter_star', a, b, format='(d,d)', /silent 
plot, a, b, /xlog, /ylog, ytitle='!6photons/Hz/Msun', xtitle='!7m!X (Hz)'

readcol, '../dndm_PopII_salpeter_Z0.004', a, b, format='(d,d)', /silent 
oplot, a, b, linestyle=2

vline, 3.288d15, linestyle=2


