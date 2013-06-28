
; File: gpi2.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.2 $) 

; Plots gamma_pi evolution along with various observational data
; points.  Explicitly compares two different Pop. III IMFs. 

;; set_plot, 'ps'
;; device, filename='gpi.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

window, xsize=1000, ysize=1000
Device, decomposed=0
!P.multi = [0,1,2,0,1]
TvLCT, 255, 0, 0, 2 
TvLCT, 255, 255, 0, 3
TvLCT, 242, 230, 56, 5
TvLCT, 127, 0, 153, 4
!P.charsize = 3
!P.thick = 1
!P.charthick = 1

; Plot our result. 
restore, 'reionfiletemplate.sav'
reiondata = read_ascii('set5/reion.out', template=reionfiletemplate)
redshift = reiondata.z
gpi = reiondata.field04
gpi = gpi * 1.0e12 
plot, redshift, gpi, /ylog, xrange=[2,100], yrange=[1.0e-14, 1.0e1], xstyle=1, xtitle='!6z', $
      ytitle='log!D10!N(!7C!6!DHI!N/10!E-12!Ns!E-1!N)',$
      ytickformat='exp2', /xlog, position=[0.1,0.1,0.9,0.7], /nodata 
gpi1 = gpi 
redshift1 = redshift 

reiondata = read_ascii('set7/reion.out', template=reionfiletemplate)
redshift = reiondata.z
gpi2 = reiondata.field04
gpi2 = gpi2 * 1.0e12 
oplot, redshift, gpi2, linestyle=5, color=2
oplot, redshift1, gpi1 

ratio = gpi/gpi2

; Plot Meiksin and White points.
plotsym, 0, 1.5, /FILL
readcol, '../data/gammapi_mw.dat', x, y, dy1, dy2, /silent 
oplot, x, y, psym=8, color=3
;; oploterror, x, y, dy1, errcolor=3, psym=3, /hibar
;; oploterror, x, y, dy2, errcolor=3, psym=3, /lobar
plot_err, x, y, dy1+dy2, color=3

x0 = x[0]
y0 = y[0] 
x1 = x[0]
y1 = 7.5e-2

; Plot Bolton and Haehnelt points.
readcol, '../data/gammapi_bh.dat', x, y, dy1, dy2, /silent 
x[0] = 0.0
oplot, x, y, psym=8, color=2
;; oploterror, x, y, dy1, errcolor=2, psym=3, /hibar
;; oploterror, x, y, dy2, errcolor=2, psym=3, /lobar
plot_err, x, y, dy=dy1, dy1=dy2, color=2
x0 = x[4]
y0 = y[4] 
x1 = x[4]
y1 = 1.0e-1 
arrow, x0, y0, x1, y1, /data, hsize=7.0, color=2 

; Plot Faucher-Giguere points. 
readcol, '../data/gammapi_cafg.dat', x, y, dy, /silent 
x[0] = 0.0 
oplot, x, y, psym=8, color=4


plot_err, x, y, dy, color=4

tick = replicate(' ',3)
plot, redshift, ratio, position=[0.1,0.7,0.9,0.97], /xlog, /ylog, xrange=[2,100], $
      xstyle=1, xtickname=tick, ytitle='ratio'

