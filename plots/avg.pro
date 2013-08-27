
; File: avg.pro
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; This is a slave code for avg_multi.pro.  See that file to understand
; what this does. 

openr, lun, 'avgset3/avg_cbyfe_1_100_lag0.dat', /get_lun
cs_dat = fltarr(2,18)
readf, lun, cs_dat
rs = cs_dat[0,*]
cbyfe_avg = cs_dat[1,*]
plotsym, 0, 0.7, /FILL
plot, rs, cbyfe_avg, yrange=[-1.2,1.0], ytitle='[C/Fe]', $
      xrange=[2,6.5], xstyle=1, ystyle=1
close, lun

for i = 0, size(rs,/n_elements)-1 do begin 
   t = interpol(tarr,zarr,rs[i])
   tnew = t + 3.0e8
   znew = interpol(zarr,tarr,tnew) 
   rs[i] = znew
endfor 
oplot, rs, cbyfe_avg, color=2

for i = 0, size(rs,/n_elements)-1 do begin 
   t = interpol(tarr,zarr,rs[i])
   tnew = t + 1.0e9
   znew = interpol(zarr,tarr,tnew) 
   rs[i] = znew
endfor 
oplot, rs, cbyfe_avg, color=3

;; openr, lun, 'avgset2/avg_cbyfe_100_260_lag0.dat', /get_lun
;; cs_dat = fltarr(2,10)
;; readf, lun, cs_dat
;; rs = cs_dat[0,*]
;; cbyfe_avg = cs_dat[1,*]
;; plotsym, 0, 0.7, /FILL
;; oplot, rs, cbyfe_avg, linestyle=2
;; close, lun

;; openr, lun, 'avgset2/avg_cbyfe_1_100_lag8.dat', /get_lun
;; cs_dat = fltarr(2,10)
;; readf, lun, cs_dat
;; rs = cs_dat[0,*]
;; cbyfe_avg = cs_dat[1,*]
;; plotsym, 0, 0.7, /FILL
;; oplot, rs, cbyfe_avg, color=2
;; close, lun

;; openr, lun, 'avgset2/avg_cbyfe_100_260_lag8.dat', /get_lun
;; cs_dat = fltarr(2,10)
;; readf, lun, cs_dat
;; rs = cs_dat[0,*]
;; cbyfe_avg = cs_dat[1,*]
;; plotsym, 0, 0.7, /FILL
;; oplot, rs, cbyfe_avg, color=2, linestyle=2
;; close, lun

;; openr, lun, 'avgset2/avg_cbyfe_1_100_lag9.dat', /get_lun
;; cs_dat = fltarr(2,10)
;; readf, lun, cs_dat
;; rs = cs_dat[0,*]
;; cbyfe_avg = cs_dat[1,*]
;; plotsym, 0, 0.7, /FILL
;; oplot, rs, cbyfe_avg, color=3
;; close, lun

;; openr, lun, 'avgset2/avg_cbyfe_100_260_lag9.dat', /get_lun
;; cs_dat = fltarr(2,10)
;; readf, lun, cs_dat
;; rs = cs_dat[0,*]
;; cbyfe_avg = cs_dat[1,*]
;; plotsym, 0, 0.7, /FILL
;; oplot, rs, cbyfe_avg, color=3, linestyle=2
;; close, lun

