
; File: fb.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; This code was used to produce a figure for the Fachbeirat report. 

set_plot, 'ps'
device, filename='fb.ps', xsize=5.0, ysize=3.1, /inches, color=1, yoffset=1.0

;window, xsize=1000, ysize=618.1
!P.multi = [0, 2, 2, 0, 1]
;Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
!P.charsize = 0.6
!P.thick = 1.0

cs4, 87, 6.0, 1.0e-2, 1.0e1 
@avg2
cs3, 87, 6.0 
@avg3

device, /close_file
set_plot, 'X'
