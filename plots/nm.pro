window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
!P.charsize = 2

restore, 'halo_template.sav'
n_data = read_ascii('set17/nofmc.out', template=stars_template)
m_data = read_ascii('set17/halos.out', template=stars_template)
row = 93
lcolumn = 150
nm = n_data.field001[lcolumn:*,row]*1.0e10 
m = m_data.field001[lcolumn:*,row]*1.0e10
z = n_data.field001[0:*]
; plot, m, nm, /xlog, /ylog, yrange=[1.0e-6,1.0e-2], xrange=[1.0e9,1.0e13], xstyle=1
plot, m, nm, /xlog, /ylog

row = 92
lcolumn = 1
nm = n_data.field001[lcolumn:*,row]*1.0e10
m = m_data.field001[lcolumn:*,row]*1.0e10
z = n_data.field001[0:*]
oplot, m, nm, linestyle=2

row = 93
lcolumn = 160
nm = n_data.field001[lcolumn:*,row]*1.0e10
m = m_data.field001[lcolumn:*,row]*1.0e10
z = n_data.field001[0:*]
oplot, m, nm, color=2, thick=3



