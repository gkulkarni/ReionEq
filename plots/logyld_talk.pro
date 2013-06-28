set_plot, 'ps'
device, filename='logyld_talk.ps', xsize=7.0, ysize=7.0, $
        /inches, color=1, /HELVETICA, yoffset=1.0
!P.font=0
!P.charsize=1.5
!P.thick = 3.0
!P.charthick = 1

TvLCT, 255, 0, 0, 2 
TvLCT, 0, 255, 0, 3
TvLCT, 0, 0, 255, 4
TvLCT, 127, 127, 0, 5
TvLCT, 127, 0, 127, 6
TvLCT, 1, 3, 64, 7

readcol, 'logyld_t.dat', c, n, o, si, fe, zn, /silent 
c = 10.0^c 
n = 10.0^n 
o = 10.0^o 
si = 10.0^si
fe = 10.0^fe
zn = 10.0^zn 

tick = [' ', '[1-100]', ' ', '[35-100]', ' ', '[100-260]', ' ']
plot, c, psym=-6, /ylog, yrange=[1.0e-5,1.0], xthick=3, ythick=3, $
      xrange=[-0.1,2.1], ytitle='yield', xtitle='IMF', xminor=1, $
      xtickname=tick, /nodata
oplot, c, psym=-6, color=2
oplot, n, psym=-6, color=3
oplot, o, psym=-6, color=4
oplot, si, psym=-6, color=5
oplot, fe, psym=-6, color=6
oplot, zn, psym=-6, color=7

xyouts, 0.5, 1.2e-2, 'C', color=2
xyouts, 0.5, 1.0e-3, 'Si', color=5
xyouts, 0.5, 8.0e-2, 'O', color=4
xyouts, 0.5, 1.0e-4, 'Fe', color=6

device, /close_file
set_plot, 'X'
