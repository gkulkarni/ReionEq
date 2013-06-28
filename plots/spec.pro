
; File: spec.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Plots the Starburst99 spectrum used in reion, just as a check. 

window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 127, 0, 153, 4
!P.charsize = 1.5
!P.thick = 1
!P.charthick = 1

speed_of_light = 2.998d10 ; cm/s
readcol, '../popsyn/spectrum', t, lambda, lum 
wl = extrac(lambda, 0, 1220) ; Ang 
l = extrac(lum, 0, 1220) ; erg/s/Ang
l = double(10.0^l) ; erg/s/Ang
wl = double(wl)
freq = dblarr(1220)
lnu = dblarr(1220) 

plot, wl, l, /xlog, /ylog
for i = 0, 1219 do begin 
   freq[i] = speed_of_light/(wl[i]*1.0e-8) ; Hz 
   lnu[i] = l[i]*wl[i]/freq[i] ; erg/s/Hz 
   
endfor
; plot, freq, lnu, /xlog, /ylog

END





