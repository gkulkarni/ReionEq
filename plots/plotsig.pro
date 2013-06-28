; File: plotsig.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Plots the D-statistic and associated p-values from 100's of runs of
; ks.py.

sig = fltarr(2,100)
openr, lun, '../foo4', /get_lun
readf, lun, sig
d = sig[0,*]
p = sig[1,*]
; plothist, p, bin=0.00001
; plothist, p, bin=0.0001, xrange=[0.0,0.005]

