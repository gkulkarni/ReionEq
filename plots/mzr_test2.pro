window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4
!P.charsize = 2

row = 80
lcolumn = 2

restore, 'halo_template.sav'
stars_data = read_ascii('set10/halos_stars.out', template=stars_template)
o_data = read_ascii('set10/halos_o.out', template=stars_template)
gas_data =  read_ascii('set10/halos_gas.out', template=stars_template)

stars = stars_data.field001[lcolumn:*,row]
gas = gas_data.field001[lcolumn:*,row]
o = o_data.field001[lcolumn:*,row]

plot, gas, /ylog

stars_data = read_ascii('set49/halos_stars.out', template=stars_template)
o_data = read_ascii('set49/halos_o.out', template=stars_template)
gas_data =  read_ascii('set49/halos_gas.out', template=stars_template)

stars = stars_data.field001[lcolumn:*,row]
gas = gas_data.field001[lcolumn:*,row]
o = o_data.field001[lcolumn:*,row]

oplot, gas, color=2







