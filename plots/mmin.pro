; File: mmin.pro 
;  Cre: 2012
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $)

; Plots the minimum mass of star forming haloes that results from
; photoionization feedback in reion. 


set_plot, 'ps'
device, filename='mmin0.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

;; window, xsize=1000, ysize=1000
;; Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
!P.charsize = 1.5
!P.thick = 1.0
!P.charthick = 1

restore, 'mmintemplate.sav'
mmin_data = read_ascii('set5/mmin.out', template=mmintemplate) 
z = mmin_data.field1
mminc = mmin_data.field2*1.0e10
mminh = mmin_data.field3*1.0e10
mminav = mmin_data.field4*1.0e10

sol = "156B
sun = '!9' + string(sol) + '!X'
plot, z, mminav, /xlog, /ylog, xrange=[2,100], xtitle='!6z', ytitle='M!Dmin!N (M!D' + sun + '!N)', xstyle=1
oplot, z, mminc, linestyle=5, color=3
oplot, z, mminh, linestyle=5, color=2

vline, 8.5, linestyle=2
xyouts, 7.5, 7.0e6, 'z!Dreion!N', orientation=90.0, charsize=2.0, alignment=0.5

legend, ['HI regions', 'HII regions', 'average'], linestyle=[5,5,0], color=[3,2,-1], /bottom, charsize=1


device, /close_file
set_plot, 'X'

