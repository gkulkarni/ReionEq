
; File: temp_talk.pro 
;  Cre: 2013-02-06
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Plot temperature evolution. 

set_plot, 'ps'
device, filename='t0_talk.ps', xsize=7.0, ysize=7.0, $
        /inches, color=1, /HELVETICA, yoffset=1.0
!P.font = 0 

TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4
!P.charsize = 1.5
!P.thick = 3
!P.charthick = 1

restore, 'reionfiletemplate.sav'
reiondata = read_ascii('set49/reion.out', template=reionfiletemplate)
redshift = reiondata.z
q = reiondata.q
tempc = reiondata.field06
temph = reiondata.field05
tempav = reiondata.field07
temphvolav = reiondata.field14
tva = temphvolav 

for i = 0, 99 do begin 
   
   tva[i] = q[i]*temphvolav[i] + (1.0-q[i])*tempc[i] 

endfor

plot, redshift, tva, /xlog, /ylog, xrange=[1,100], xtitle='redshift', $
      ytitle='T!D0!N (K)', ythick=3, xthick=3 
oplot, redshift, tempc, linestyle=2, color=1, thick=3 
oplot, redshift, temphvolav, linestyle=2, color=2, thick=3

legend, ['H I', 'H II'], linestyle=[2,2], color=[3,2], $
        /bottom, box=0

;; vline, 7.0, linestyle=2
;; xyouts, 6.2, 1.0e3, 'z!Dreion!N', orientation=90.0, charsize=2.0, alignment=0.5

device, /close_file
set_plot, 'X'

END
