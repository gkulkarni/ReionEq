
; File: samples.pro 
;  Cre: 01-2013 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Plots the random sample produced by sample_ns.f90.  Also plots the
; DLA distribution in an adjacent panel. 

set_plot, 'ps'
device, filename='sample.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

;; window, xsize=1000, ysize=1000
;; Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 127, 0, 153, 4
!P.charsize = 1.5
!P.thick = 1.0
; !P.multi=[1,1,2,0,0]

array_size = 10 

erase 
multiplot, [1,2], mXtitle='!6[O/Si]', mXtitsize='1.5', $
           mytitsize='1.5', mYtitoffset=1.5, mXtitoffset=1.5, $
           /doyaxis
ns_p2 = fltarr(2,176)
openr, lun, 'ns_p2.dat', /get_lun
readf, lun, ns_p2
obysi_p2 = ns_p2[0,*] 
ndla_p2 = ns_p2[1,*]
plot, obysi_p2, ndla_p2, ytitle='d!E2!NN / dXd[O/Si]', $
      xrange=[0.0,0.6]
close, lun

ns_p3 = fltarr(2,88)
openr, lun, 'ns_p3.dat', /get_lun
readf, lun, ns_p3
obysi_p3 = ns_p3[0,*] 
ndla_p3 = ns_p3[1,*]
oplot, obysi_p3, ndla_p3, color=2
close, lun

;' legend, [label, 'Fontanot et al. 2007'], linestyle=[-1,0], /right, box=0

legend, ['z=6'], /right, box=0, charsize=2

multiplot 

np3_sample = fltarr(array_size)
openr, lun, '../noisy_sample_p3.dat', /get_lun
readf, lun, np3_sample
plothist, np3_sample, bin=0.05, ytitle='number of DLAs', $
          xrange=[0.0,0.6], color=2, yrange=[0,4]
close, lun
multiplot, /reset 
close, /all

legend, ['z=6'], /right, box=0, charsize=2

device, /close_file
set_plot, 'X'

