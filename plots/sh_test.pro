window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4
!P.charsize = 1.5
!P.thick = 1.0

erase
multiplot, [1,2], mXtitle='!6z', mXtitsize='1.5', mytitsize='1.5', $
           mYtitoffset=1.5, mXtitoffset=1.5, /doyaxis

restore, 'mmintemplate.sav'
mmin_data = read_ascii('set7/mmin.out', template=mmintemplate) 
z = mmin_data.field1
mminc = mmin_data.field2
mminh = mmin_data.field3
mminav = mmin_data.field4

mminh=alog10(mminh)+10.0
mminav=alog10(mminav)+10.0
plot, z, mminh, /xlog, xrange=[1,100], linestyle=5, ytitle = 'log!D10!N(M/M!D!9n!X!N)', $
      yrange=[6,12]
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

xyouts, 30.0, 10.0, 'model 3', charsize=2.0, alignment=0.5

multiplot 

;---------------------------------------

solar = 1.17
o_data = read_ascii('set7/halos_o.out', template=stars_template)
si_data =  read_ascii('set7/halos_Si.out', template=stars_template)

o = o_data.field001[100,*]
si = si_data.field001[100,*]
obysi = fltarr(100)
for i = 0, 99 do begin 
   if (o[0,i] eq 0.0) then begin 
      obysi[i] = -10.0
   endif else begin 
      obysi[i] = alog10(abs(o[0,i]/si[0,i])) - solar 
   endelse
endfor 
plot, z, obysi, /xlog, xrange=[1,100], yrange=[-1,1], ytitle='[O/Si]', xstyle=1
; oplot, z, obysi, color=2

o = o_data.field001[140,*]
si = si_data.field001[140,*]
obysi = fltarr(100)
for i = 0, 99 do begin 
   if (o[0,i] eq 0.0) then begin 
      obysi[i] = -10.0
   endif else begin 
      obysi[i] = alog10(abs(o[0,i]/si[0,i])) - solar 
   endelse
endfor 
oplot, z, obysi, color=3

o = o_data.field001[170,*]
si = si_data.field001[170,*]
obysi = fltarr(100)
for i = 0, 99 do begin 
   if (o[0,i] eq 0.0) then begin 
      obysi[i] = -10.0
   endif else begin 
      obysi[i] = alog10(abs(o[0,i]/si[0,i])) - solar 
   endelse
endfor 
oplot, z, obysi, color=2
vline, 6.0, linestyle=2

multiplot, /reset 
close, /all 

END
