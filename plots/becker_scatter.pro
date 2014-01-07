set_plot, 'ps'
device, filename='scatter.ps', xsize=7.0, ysize=7.0, /inches, $
        color=1, yoffset=1.0

TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
!P.charsize = 1.5
!P.thick = 1.0
!P.charthick = 1

; Read [C/Fe].
readcol, 'cbyfe.dat', z, cbyfe, cbyfe_err, format='(F,F,F)'

; Read [O/Si]. 
readcol, 'obysi.dat', z, obysi, obysi_err, format='(F,F,F)'
obysi = -obysi 

; Do the plotting.
plot, cbyfe, obysi, psym=6, xrange=[-0.5,0.5], yrange=[-0.5,0.5], $
      xtitle='!6[C/Fe]', ytitle='[O/Si]'

device, /close_file
set_plot, 'X'

END

