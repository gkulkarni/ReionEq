row = 99

stars = stars_data.field001[lcolumn:*,row]
gas = gas_data.field001[lcolumn:*,row]
o = o_data.field001[lcolumn:*,row]

array_size = size(stars)
array_length = array_size[1] 
met = fltarr(array_length)
obyh = fltarr(array_length)

for i = 0, array_length-1 do begin 

   mH = gas[i]*0.71

   if (o[i] eq 0.0) then begin 
      obyh[i] = 0.0
   endif else begin 
      obyh[i] = alog10(abs(o[i]/mH)) + 1.87
   endelse

   met[i] = obyh[i] - 1.87 + 12.0 

endfor

mets = smooth(met, 30, /edge_truncate) 
stars = alog10(stars)  + 10.0

oplot, stars, mets, color=2

readcol, '../data/mzr_tremonti.dat', mst, metlow, met, methigh 
metlow = met - metlow
methigh = methigh - met 
oplot, mst, met, psym=8, color=2
plot_err, mst, met, metlow, dy2=methigh, color=2

;-------------------------------

row = 92

stars = stars_data.field001[lcolumn:*,row]
gas = gas_data.field001[lcolumn:*,row]
o = o_data.field001[lcolumn:*,row]

array_size = size(stars)
array_length = array_size[1] 
met = fltarr(array_length)
obyh = fltarr(array_length)

for i = 0, array_length-1 do begin 

   mH = gas[i]*0.71

   if (o[i] eq 0.0) then begin 
      obyh[i] = 0.0
   endif else begin 
      obyh[i] = alog10(abs(o[i]/mH)) + 1.87
   endelse

   met[i] = obyh[i] - 1.87 + 12.0 - 1.1

endfor

mets = smooth(met, 30, /edge_truncate) 
stars = alog10(stars) + 10.0 - 1.0

oplot, stars, mets, color=2

readcol, '../data/mzr_maiolino.dat', mst, dx2, dx1, met, dy2, dy1
oplot, mst, met, psym=8, color=2
plot_err, mst, met, dy1, dy2=dy2, dx1=dx1, dx2=dx2, color=2


3.56270e+07      50.0093
5.97869e+09      1.00003
