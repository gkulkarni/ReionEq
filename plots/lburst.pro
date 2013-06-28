window, xsize=1000, ysize=1000
!p.charsize=2
readcol, '../popsyn/l1500-spectrum-lt.dat', age, lamb, lum
plot, age, lum, /xlog
