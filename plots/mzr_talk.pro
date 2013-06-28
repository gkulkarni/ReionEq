
; File: mzr_talk.pro
;  Cre: 2013-02-06
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.3 $)

set_plot, 'ps'
device, filename='mzr_talk.ps', xsize=7.0, ysize=7.0, $
        /inches, color=1, /HELVETICA, yoffset=1.0
!P.font = 0 

TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4
!P.charsize = 1.5
!P.thick = 3
!P.charthick = 1

row = 94
lcolumn = 95

restore, 'halo_template.sav'
stars_data = read_ascii('set49/halos_stars.out', template=stars_template)
o_data = read_ascii('set49/halos_o.out', template=stars_template)
gas_data =  read_ascii('set49/halos_gas.out', template=stars_template)

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

   met[i] = obyh[i] - 1.87 + 12.0 - 1.0

endfor

mets = smooth(met, 30, /edge_truncate) 
stars = alog10(stars)  + 10.0 - 1.0

plot, stars, mets, yrange=[7,9], xtitle='log!D10!N(M!D*!N/Msun)', $
      ytitle='12+log!D10!N(O/H)', xrange=[8,12], xthick=3, ythick=3

readcol, '../data/erb_pointscorrect.dat', str, met 
plotsym, 0, 1, /FILL
oplot, str, met, psym=8 

readcol, '../data/erb_xlow.dat', strlow, foo
readcol, '../data/erb_xhigh.dat', strhigh, foo
readcol, '../data/erb_ylow.dat', foo, metlow 
readcol, '../data/erb_yhigh.dat', foo, methigh

strlow = str - strlow 
strhigh = strhigh - str 
metlow = met - metlow 
methigh = methigh - met 

plot_err, str, met, metlow, dy2=methigh, dx1=strlow, dx2=strhigh

;-------------------------------

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

; oplot, stars, mets, color=2

readcol, '../data/mzr_tremonti.dat', mst, metlow, met, methigh 
metlow = met - metlow
methigh = methigh - met 
; oplot, mst, met, psym=8, color=2
; plot_err, mst, met, metlow, dy2=methigh, color=2

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

legend, ['z = 2.3', 'z = 3.7'], linestyle=[0,0], color=[1,2], $
        /bottom, /right, box=0

device, /close_file
set_plot, 'X'

END

