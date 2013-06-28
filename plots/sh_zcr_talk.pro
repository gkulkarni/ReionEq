
; File: sh_zcr.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Plots mass and metallicity evolution of individual haloes from
; reion.  This plot was used in paper. 

;; set_plot, 'ps'
;; device, filename='sh_zcr.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

set_plot, 'ps'
device, filename='sh_zcr_talk.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

!P.font = 0 
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4
!P.charsize = 1.5
!P.thick = 3
!P.charthick=1

erase
multiplot, [1,2], mXtitle='z', mXtitsize='1.5', mytitsize='1.5', $
           mYtitoffset=1.5, mXtitoffset=1.5, /doyaxis

restore, 'mmintemplate.sav'
mmin_data = read_ascii('set5/mmin.out', template=mmintemplate) 
z = mmin_data.field1
mminc = mmin_data.field2
mminh = mmin_data.field3
mminav = mmin_data.field4
mminh=alog10(mminh)+10.0
mminav=alog10(mminav)+10.0
plot, z, mminh, /xlog, xrange=[1,100], linestyle=5, ytitle = 'log!D10!N(M / Msun)', yrange=[6,12]

restore, 'halo_template.sav'
halos_data = read_ascii('set5/halos.out', template=stars_template)
halos = halos_data.field001[170,*]
halos=alog10(halos)+10.0
oplot, z, halos, color=2
halos = alog10(halos_data.field001[140,*])+10.0
oplot, z, halos, color=3
halos = alog10(halos_data.field001[100,*])+10.0
oplot, z, halos

; xyouts, 30.0, 10.0, 'model 3', charsize=2.0, alignment=0.5
;; legend, ['Pop. III IMF: 100-260 M!D'+sun], box=0, position=[9.7,11.5]
;; legend, ['!NM!Dmin!N (HII region)'], linestyle=5, position=[8.5,11.0], box=0
 multiplot 

;---------------------------------------

metal_data = read_ascii('set5/halos_metals.out', template=stars_template)
gas_data =  read_ascii('set5/halos_gas.out', template=stars_template)
metal = metal_data.field001[100,*]
gas = gas_data.field001[100,*]
Zmetal = fltarr(100)
for i = 0, 99 do begin 
   Zmetal[i] = (metal[0,i]/gas[0,i])/0.02
endfor 
; plot, z, Zmetal, /xlog, xrange=[1,100], yrange=[-1,1], ytitle='Z',
; xstyle=1
plot, z, Zmetal, /xlog, xrange=[1,100], ytitle='Z / Zsun', $
      xstyle=1, /ylog, ytickformat='Exp1', yrange=[1.0e-6,1.0]

metal_data = read_ascii('set5/halos_metals.out', template=stars_template)
gas_data =  read_ascii('set5/halos_gas.out', template=stars_template)
metal = metal_data.field001[140,*]
gas = gas_data.field001[140,*]
Zmetal = fltarr(100)
for i = 0, 99 do begin 
   Zmetal[i] = (metal[0,i]/gas[0,i])/0.02
endfor 
oplot, z, Zmetal, color=3

metal_data = read_ascii('set5/halos_metals.out', template=stars_template)
gas_data =  read_ascii('set5/halos_gas.out', template=stars_template)
metal = metal_data.field001[170,*]
gas = gas_data.field001[170,*]
Zmetal = fltarr(100)
for i = 0, 99 do begin 
   Zmetal[i] = (metal[0,i]/gas[0,i])/0.02
endfor 
oplot, z, Zmetal, color=2

hline, 1.0e-4, linestyle=2
xyouts, 2.0, 2.0e-4, 'Z!Dcrit!N', charsize=1.5, alignment=0.5

multiplot, /reset 
close, /all 

device, /close_file
set_plot, 'X'

END
