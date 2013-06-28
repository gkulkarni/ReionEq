
window, xsize=600, ysize=600
Device, decomposed=0
!p.charsize=2

readcol, 'logyld_t.dat', c, n, o, si, fe, zn, /silent 
c = 10.0^c 
n = 10.0^n 
o = 10.0^o 
si = 10.0^si
fe = 10.0^fe
zn = 10.0^zn 

plot, c, psym=-6, /ylog, yrange=[1.0e-5,1.0-1]
oplot, n, psym=-6
oplot, o, psym=-6
oplot, si, psym=-6
oplot, fe, psym=-6
oplot, zn, psym=-6



