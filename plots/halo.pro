PRO halo, column 

set_plot, 'ps'
device, filename='halo.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

;window, xsize=500, ysize=500
;!P.multi = 0
;Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4
!P.charsize = 2.0
!P.thick = 1.0
!P.charthick = 3 
!P.thick = 3 

restore, 'halo_template.sav'
halos_data = read_ascii('set49/halos.out', template=stars_template)
z = halos_data.field001[0,*]
halos = halos_data.field001[column,*]
halos = halos*1.0e10 
;halos=alog10(halos)
plot, z, halos, /xlog, xrange=[1,100], xtitle='!6z', ytitle='mass', /ylog

stars_data = read_ascii('set49/halos_stars.out', template=stars_template)
stars = stars_data.field001[column,*]
stars = stars*1.0e10
;stars = alog10(stars) 
oplot, z, stars, color=2

gas_data = read_ascii('set49/halos_gas.out', template=stars_template)
gas = gas_data.field001[column,*]
gas = gas*1.0e10
;gas = alog10(gas) 
oplot, z, gas, color=3

metals_data = read_ascii('set49/halos_metals.out', template=stars_template)
metals = metals_data.field001[column,*]
metals = metals*1.0e10 
;metals = alog10(metals) 
oplot, z, metals, color=3, linestyle=5

device, /close_file
set_plot, 'X'

END
