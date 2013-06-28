
set_plot, 'ps'
device, filename='sed_talk.ps', xsize=7.0, ysize=7.0, $
        /inches, color=1, /HELVETICA, yoffset=1.0
!P.font = 0 

TvLCT, 200, 200, 200, 2 
!P.charsize = 1.5
!P.thick = 3
!P.charthick = 1

readcol, '../dndm_PopIII_salpeter_star', a, b, format='(d,d)', /silent 
cgplot, a, b, /xlog, /ylog, ytitle='photons / Hz / Msun', xtitle='$\nu$ (Hz)', thick=5, xrange=[2.0e15,1.0d17], xstyle=1

readcol, '../dndm_PopII_salpeter_Z0.004', a, b, format='(d,d)', /silent 
oplot, a, b, color=2, thick=5

vline, 3.288d15, linestyle=2

legend, ['Pop II', 'Pop III'], linestyle=[0,0], color=[2,1], thick=[5,5], /right, box=0

device, /close_file
set_plot, 'X'
