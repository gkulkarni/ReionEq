
; File: gpi.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.2 $) 

; Plots gamma_pi evolution along with various observational data
; points.  This figure was used in the paper. 

;; set_plot, 'ps'
;; device, filename='gpi.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
;TvLCT, 255, 255, 0, 4
TvLCT, 127, 0, 153, 4
TvLCT, 255, 255, 0, 5
TvLCT, 1, 3, 64, 6
!P.charsize = 1.5
!P.thick = 1
!P.charthick = 1
!P.multi=1

restore, 'reionfiletemplate.sav'
reiondata = read_ascii('set48/reion.out', template=reionfiletemplate)
redshift = reiondata.z
gpi = reiondata.field04
gpi = gpi * 1.0e12 * 0.1
ticks = [2,3,4,5,6,7,8,9,10,20]
nticks = n_elements(ticks)
plot, redshift, gpi, /ylog, xrange=[2,30], xstyle=1, xtitle='!6z', $
      ytitle='log!D10!N(!7C!6!DHI!N/10!E-12!Ns!E-1!N)', $
      yrange=[1.0e-4, 1.0e1], ytickformat='exponent', /xlog , $
      xticks=nticks-1, xtickv=ticks

reiondata = read_ascii('set27/reion.out', template=reionfiletemplate)
redshift = reiondata.z
gpi = reiondata.field04
gpi = gpi * 1.0e12 
oplot, redshift, gpi, linestyle=5

openr, lun, '../bh.dat', /get_lun
obsdata = fltarr(4,3)
readf, lun, obsdata
zobs = obsdata[0,*]
gobs = obsdata[1,*]*1.0e12 
gerr_low = obsdata[2,*]*1.0e12 
gerr_high = obsdata[3,*]*1.0e12 
gerr_low[2] = 5.0e-2
gerr_low = gobs - gerr_low 
gerr_high[2] = gerr_high[2]+1.0e-13
gerr_high = gerr_high - gobs

plotsym, 0, 1, /FILL
zerr = fltarr(3)

restore, 'reionfiletemplate.sav'
reiondata = read_ascii('set26/reion.out', template=reionfiletemplate)
redshift = reiondata.z
gpi = reiondata.field04

readcol, '../data/gammapi_mw.dat', x, y, dy1, dy2 
oplot, x, y, psym=8, color=3
oploterror, x, y, dy1, errcolor=3, psym=3, /hibar
oploterror, x, y, dy2, errcolor=3, psym=3, /lobar
x0 = x[0]
y0 = y[0] 
x1 = x[0]
y1 = 7.5e-2
;; arrow, x0, y0, x1, y1, /data, hsize=7.0, color=3

readcol, '../data/gammapi_bh.dat', x, y, dy1, dy2 
x[0] = 0.0
oplot, x, y, psym=8, color=2
oploterror, x, y, dy1, errcolor=2, psym=3, /hibar
oploterror, x, y, dy2, errcolor=2, psym=3, /lobar

x0 = x[4]
y0 = y[4] 
x1 = x[4]
y1 = 1.0e-1 
arrow, x0, y0, x1, y1, /data, hsize=7.0, color=2 

readcol, '../data/gammapi_cafg.dat', x, y, dy
x[0] = 0.0 
oplot, x, y, psym=8, color=4
oploterror, x, y, dy, errcolor=4, psym=3

vline, 8.5, linestyle=2
xyouts, 8.0, 1.0e0, 'z!Dreion!N', orientation=90.0, charsize=2.0, alignment=0.5

readcol, 'gpi_old5.dat', z, q, temph, dnlldz, gammapi, nphdot, lmfp, source
; plot, z, gammapi/1.0e-12, /xlog, /ylog
; oplot, z, gammapi/1.0e-12, color=5, linestyle=2
;oplot, z, gammapi/1.0e-12, color=6, linestyle=2

; Plot Tirth result. 
readcol, 'WMAP5atomic_IGM.dat', z, q, gpi1, gpi2, gpi3, t1, t2, $
         t3, t4, dn, qh, xhi, ng, t5
;oplot, z, gpi1/1.0e-12 

legend, ['Faucher-Giguere 08', 'Meiksin and White 04', 'Bolton and Haehnelt 07'], linestyle=[0,0,0], color=[4,3,2], psym=[8,8,8], /bottom, charsize=1.2

;; device, /close_file
;; set_plot, 'X'


