
; File: avg.pro
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; This is a slave code for avg_multi.pro.  See that file to understand
; what this does. 

openr, lun, 'avgset4/avg_cbyfe_1_100_f.dat', /get_lun
cs_dat = fltarr(2,18)
readf, lun, cs_dat
rs = cs_dat[0,*]
cbyfe_avg = cs_dat[1,*]
plotsym, 0, 0.7, /FILL
plot, rs, cbyfe_avg, yrange=[-1.2,1.0], ytitle='[C/Fe]', $
      xrange=[2,10.5], xstyle=1, ystyle=1
close, lun

openr, lun, 'avgset4/avg_cbyfe_1_100_m.dat', /get_lun
cs_dat = fltarr(2,18)
readf, lun, cs_dat
rs = cs_dat[0,*]
cbyfe_avg = cs_dat[1,*]
plotsym, 0, 0.7, /FILL
oplot, rs, cbyfe_avg, color=2
close, lun

openr, lun, 'avgset4/avg_cbyfe_1_100_o.dat', /get_lun
cs_dat = fltarr(2,18)
readf, lun, cs_dat
rs = cs_dat[0,*]
cbyfe_avg = cs_dat[1,*]
plotsym, 0, 0.7, /FILL
oplot, rs, cbyfe_avg, color=3
close, lun

