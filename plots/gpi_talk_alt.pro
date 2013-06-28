
set_plot, 'ps'
device, filename='gpi_talk_alt.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

!P.font = 0 
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 242, 230, 56, 5
TvLCT, 127, 0, 153, 4
!P.charsize = 1.5
!P.thick = 3
!P.charthick = 1

; Plot our result. 
restore, 'reionfiletemplate.sav'
reiondata = read_ascii('set5/reion.out', template=reionfiletemplate)
redshift = reiondata.z
gpi = reiondata.field04
gpi = gpi * 1.0e12 
cgplot, redshift, gpi, /ylog, xrange=[2,10], yrange=[0.001,2.0], xstyle=1, xtitle='z', $
        ytitle='log!D10!N($\Gamma$!DHI!N / 10!E-12!Ns!E-1!D)', /xlog, /nodata 
gpi1 = gpi 
redshift1 = redshift 

reiondata = read_ascii('set7/reion.out', template=reionfiletemplate)
redshift = reiondata.z
gpi2 = reiondata.field04
gpi2 = gpi2 * 1.0e12 
oplot, redshift, gpi2, linestyle=5, color=2, thick=5
oplot, redshift1, gpi1, thick=5 
ratio = gpi/gpi2

; Plot Meiksin and White points.
plotsym, 0, 0.8, /FILL
readcol, '../data/gammapi_mw.dat', x, y, dy1, dy2, /silent 
oplot, x, y, psym=8, color=3
oploterror, x, y, dy1, errcolor=3, psym=3, /hibar
oploterror, x, y, dy2, errcolor=3, psym=3, /lobar
x0 = x[0]
y0 = y[0] 
x1 = x[0]
y1 = 7.5e-2

; Plot Bolton and Haehnelt points.
readcol, '../data/gammapi_bh.dat', x, y, dy1, dy2, /silent 
x[0] = 0.0
oplot, x, y, psym=8, color=2
oploterror, x, y, dy1, errcolor=2, psym=3, /hibar
oploterror, x, y, dy2, errcolor=2, psym=3, /lobar
x0 = x[4]
y0 = y[4] 
x1 = x[4]
y1 = 1.0e-1 
arrow, x0, y0, x1, y1, /data, hsize=7.0, color=2 

; Plot Faucher-Giguere points. 
readcol, '../data/gammapi_cafg.dat', x, y, dy, /silent 
x[0] = 0.0 
oplot, x, y, psym=8, color=4
oploterror, x, y, dy, errcolor=4, psym=3

legend, ['Faucher-Giguere 08', 'Meiksin and White 04', 'Bolton and Haehnelt 07'], $
        linestyle=[0,0,0], psym=[8,8,8], color=[4,3,2], /bottom, box=0

vline, 8.5, linestyle=2

device, /close_file
set_plot, 'X'
