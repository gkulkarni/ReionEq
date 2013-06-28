restore, 'mmintemplate.sav'
mmin_data = read_ascii('set69/mmin.out', template=mmintemplate) 
z = mmin_data.field1
mminc = mmin_data.field2
mminh = mmin_data.field3
mminav = mmin_data.field4
plot, z, mminh, /xlog, /ylog, xrange=[1,100] 

mmin_data = read_ascii('set67/mmin.out', template=mmintemplate) 
z = mmin_data.field1
mminc = mmin_data.field2
mminh = mmin_data.field3
mminav = mmin_data.field4
oplot, z, mminh, color=2 

mmin_data = read_ascii('set68/mmin.out', template=mmintemplate) 
z = mmin_data.field1
mminc = mmin_data.field2
mminh = mmin_data.field3
mminav = mmin_data.field4
oplot, z, mminh, color=3 

