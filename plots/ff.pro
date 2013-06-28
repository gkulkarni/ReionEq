
; File: ff.pro 
;  Cre: 2012
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $)

; Plot the HII filling factor evolution. 

;; set_plot, 'ps'
;; device, filename='ff.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 128, 128, 4
TvLCT, 0, 255, 0, 5
TvLCT, 255, 128, 0, 6
TvLCT, 128, 0, 128, 7
TvLCT, 255, 0, 255, 8
TvLCT, 255, 255, 0, 9
!P.charsize = 1.5

restore, 'reionfiletemplate.sav'
; reiondata = read_ascii('set89/reion.out',
; template=reionfiletemplate)
; reiondata = read_ascii('set5/reion.out', template=reionfiletemplate)
reiondata = read_ascii('set64/reion.out', template=reionfiletemplate)
redshift = reiondata.z
ff = reiondata.q
ff_h1 = 1.0-ff 
plot, redshift, ff_h1, xtitle='!6z', ytitle='1-Q!DHII!N', xrange=[4.5,12]

plotsym, 0, 1, /FILL
readcol, '../data/robertson2.dat', x, y, dym, dyp, ci, format='(f,f,f,f,i)', /silent 

; oplot, x, y, psym=8, color=4
n_data = size(x, /n_elements)

for i = 0, n_data-1 do begin 
   if (dym(i) eq 0.0) then begin  
      xl = fltarr(1) 
      yl = fltarr(1) 
      xl(0) = x(i) 
      yl(0) = y(i) 
      oplot, xl, yl, psym=8, color=ci(i)
      x0 = x(i) 
      y0 = y(i) 
      x1 = x(i) 
      y1 = y(i) + dyp(i) 
      arrow, x0, y0, x1, y1, /data, hsize=20.0, color=ci(i)
   endif else if (dyp(i) eq 0.0) then begin 
      xl = fltarr(1) 
      yl = fltarr(1) 
      xl(0) = x(i) 
      yl(0) = y(i) 
      oplot, xl, yl, psym=8, color=ci(i)
      x0 = x(i)
      y0 = y(i)
      x1 = x(i)
      y1 = y(i) + dym(i) 
      arrow, x0, y0, x1, y1, /data, hsize=20.0, color=ci(i)
   endif
endfor 

;; xl = fltarr(1) 
;; yl = fltarr(1) 
;; dy = yl 
;; dy2 = yl 

;; xl(0) = x(n_data-1) 
;; yl(0) = y(n_data-1) 
;; dy(0) = -dym(n_data-1)
;; dy2(0) = dyp(n_data-1)
;; plot_err, xl, yl, dy, dy2=dy2, color=4 

legend, ['Quasar near zone', 'Dark Ly!7a!X forest pixels', 'GRB damping wing',  $
         'Ly!7a!X galaxy clustering', 'Ly!7a!X emitters'], $
        linestyle=[0,0,0,0,0], psym=[8,8,8,8,8], color=[2,3,4,8,5], charsize=1.2


;; device, /close_file
;; set_plot, 'X'

END

