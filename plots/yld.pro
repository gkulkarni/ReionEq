
; File: yld.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Plots yields from yield tables. Just a visualisation aid. 

PRO yld, column, row 

; column: column number in, e.g., yields.wwhw.Z0 (begin with 1) 
;    row: row number in , e.g., yields.wwhw.Z0 (begin with 0) 

colstr = string(column, format='(I2.2)')

field_0 = 'data_0.field' + colstr
field_01 = 'data_01.field' + colstr
field_001 = 'data_001.field' + colstr
field_00001 = 'data_00001.field' + colstr

window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4
!P.charsize = 2.0
!P.thick = 1.0

restore, 'ztemplate.sav'
data_0 = read_ascii('../data/yields.wwhw.Z0', template=ztemplate)
data_01 = read_ascii('../data/yields.wwhw.Z01', template=ztemplate)
data_001 = read_ascii('../data/yields.wwhw.Z001', template=ztemplate)
data_00001 = read_ascii('../data/yields.wwhw.Z00001', template=ztemplate)

yldmass = fltarr(4)

yldmass[0] = data_0.field01[8]
yldmass[1] = data_01.field01[8]
yldmass[2] = data_001.field01[8]
yldmass[3] = data_00001.field01[8]

;; yldmass[0] = field_0[13]
;; yldmass[1] = field_01[13]
;; yldmass[2] = field_001[13]
;; yldmass[3] = field_00001[13]

zarr = fltarr(4)

zarr[0] = 0.0
zarr[1] = 0.1 
zarr[2] = 0.01
zarr[3] = 0.0001

plotsym, 0, 2, /FILL

plot, zarr, yldmass, psym=8, xrange=[-0.01, 0.12], xtitle='!6Z', ytitle='Element mass (Msun)', /ylog 
oplot, [zarr[0]], [yldmass[0]], color=2, psym=8

END
