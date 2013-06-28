
window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4
!P.charsize = 1.5
!P.thick = 1.0
!P.charthick = 1

row = 70
lcolumn = 65

restore, 'halo_template.sav'
stars_data = read_ascii('set21/halos_stars.out', template=stars_template)
o_data = read_ascii('set21/halos_o.out', template=stars_template)
gas_data =  read_ascii('set21/halos_coolgas.out', template=stars_template)

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

   met[i] = obyh[i] + 12.0 

endfor

; mets = smooth(met, 30, /edge_truncate) 
stars = alog10(stars) + 10.0

plot, stars, met, yrange=[8,14], xtitle='!6log!D10!N(M!D*!N/M!D!9n!X!N)', $
      ytitle='!612+log!D10!N(O/H)', xrange=[7,14], /nodata 

for row = 92, 99 do begin 

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

      met[i] = obyh[i] + 12.0 

   endfor

   ; mets = smooth(met, 30, /edge_truncate) 
   stars = alog10(stars) + 10.0

   if row eq 96 then begin 
      oplot, stars, met, color=2
   endif else begin 
      oplot, stars, met
   endelse

endfor 

END

