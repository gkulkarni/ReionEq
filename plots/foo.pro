restore, 'mmintemplate.sav'
mmin_data = read_ascii('set7/mmin.out', template=mmintemplate) 
z = mmin_data.field1
mminc = mmin_data.field2
mminh = mmin_data.field3
mminav = mmin_data.field4

mminh=alog10(mminh)+10.0
mminav=alog10(mminav)+10.0
plot, z, mminh, /xlog, xrange=[1,100], linestyle=5, ytitle = 'log!D10!N(M/M!D!9n!X!N)', yrange=[6,12]
oplot, z, mminav, linestyle=1

restore, 'halo_template.sav'
halos_data = read_ascii('set7/halos.out', template=stars_template)

halos = halos_data.field001[170,*]
halos=alog10(halos)+10.0
oplot, z, halos, color=2

halos = alog10(halos_data.field001[140,*])+10.0
oplot, z, halos, color=3

halos = alog10(halos_data.field001[100,*])+10.0
oplot, z, halos

