
; File: sh_obysi.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Plot evolution of [O/Si] for three different haloes in all three models.

;; set_plot, 'ps'
;; device, filename='sh_obysi.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4
TvLCT, 177, 218, 226, 5 
!P.charsize = 1.5
!P.thick = 1.0

erase
multiplot, [1,3], mXtitle='!6z', mXtitsize='1.5', mytitsize='1.5', $
           mYtitoffset=1.5, mXtitoffset=1.5, /doyaxis

halo_column = 100
restore, 'halo_template.sav'
o_data = read_ascii('set69/halos_o.out', template=stars_template)
fe_data =  read_ascii('set69/halos_Si.out', template=stars_template)
o = o_data.field001[halo_column,*]
fe = fe_data.field001[halo_column,*]
z = o_data.field001[0,*]
obyfe = fltarr(100)
for i = 0, 99 do begin 
   if (o[0,i] eq 0.0) then begin 
      obyfe[i] = -10.0
   endif else begin 
      obyfe[i] = alog10(abs(o[0,i]/fe[0,i])) - 1.17 
   endelse
endfor 
plot, z, obyfe, /xlog, xrange=[1,100], yrange=[-1,1], ytitle='[O/Si]'

o_data = read_ascii('set67/halos_o.out', template=stars_template)
fe_data =  read_ascii('set67/halos_Si.out', template=stars_template)
o = o_data.field001[halo_column,*]
fe = fe_data.field001[halo_column,*]
z = o_data.field001[0,*]
obyfe = fltarr(100)
for i = 0, 99 do begin 
   if (o[0,i] eq 0.0) then begin 
      obyfe[i] = -10.0
   endif else begin 
      obyfe[i] = alog10(abs(o[0,i]/fe[0,i])) - 1.17 
   endelse
endfor 
oplot, z, obyfe, linestyle = 3

o_data = read_ascii('set68/halos_o.out', template=stars_template)
fe_data =  read_ascii('set68/halos_Si.out', template=stars_template)
o = o_data.field001[halo_column,*]
fe = fe_data.field001[halo_column,*]
z = o_data.field001[0,*]
obyfe = fltarr(100)
for i = 0, 99 do begin 
   if (o[0,i] eq 0.0) then begin 
      obyfe[i] = -10.0
   endif else begin 
      obyfe[i] = alog10(abs(o[0,i]/fe[0,i])) - 1.17 
   endelse
   print, i, o[0,i], fe[0,i], obyfe[i]
endfor 
oplot, z, obyfe, linestyle = 1

sol = "156B
sun = '!9' + string(sol) + '!X'
xyouts, 1.2, -0.7, 'M!Dhalo!N(z=0) = 10!E9.5!NM!D' + sun, charsize=1.5

x = fltarr(2) 
y = fltarr(2) 
x(0) = 15.0
y(0) = 0.5
x(1) = 30.0
y(1) = 0.5
oplot, x, y, linestyle=1

x(0) = 15.0
y(0) = 0.3
x(1) = 30.0
y(1) = 0.3
oplot, x, y, linestyle=3

x(0) = 15.0
y(0) = 0.1
x(1) = 30.0
y(1) = 0.1
oplot, x, y

xyouts, 35.0, 0.5, 'model 1'
xyouts, 35.0, 0.3, 'model 2'
xyouts, 35.0, 0.1, 'model 3'
xyouts, 25.0, 0.7, '(all panels)', alignment=0.5



multiplot
;----------------------------------------------------

halo_column = 140
restore, 'halo_template.sav'
o_data = read_ascii('set69/halos_o.out', template=stars_template)
fe_data =  read_ascii('set69/halos_Si.out', template=stars_template)
o = o_data.field001[halo_column,*]
fe = fe_data.field001[halo_column,*]
z = o_data.field001[0,*]
obyfe = fltarr(100)
for i = 0, 99 do begin 
   if (o[0,i] eq 0.0) then begin 
      obyfe[i] = -10.0
   endif else begin 
      obyfe[i] = alog10(abs(o[0,i]/fe[0,i])) - 1.17 
   endelse
endfor 
plot, z, obyfe, /xlog, xrange=[1,100], yrange=[-1,1], ytitle='[O/Si]'
oplot, z, obyfe, color=3

o_data = read_ascii('set67/halos_o.out', template=stars_template)
fe_data =  read_ascii('set67/halos_Si.out', template=stars_template)
o = o_data.field001[halo_column,*]
fe = fe_data.field001[halo_column,*]
z = o_data.field001[0,*]
obyfe = fltarr(100)
for i = 0, 99 do begin 
   if (o[0,i] eq 0.0) then begin 
      obyfe[i] = -10.0
   endif else begin 
      obyfe[i] = alog10(abs(o[0,i]/fe[0,i])) - 1.17 
   endelse
endfor 
oplot, z, obyfe, color=3, linestyle=3

o_data = read_ascii('set68/halos_o.out', template=stars_template)
fe_data =  read_ascii('set68/halos_Si.out', template=stars_template)
o = o_data.field001[halo_column,*]
fe = fe_data.field001[halo_column,*]
z = o_data.field001[0,*]
obyfe = fltarr(100)
for i = 0, 99 do begin 
   if (o[0,i] eq 0.0) then begin 
      obyfe[i] = -10.0
   endif else begin 
      obyfe[i] = alog10(abs(o[0,i]/fe[0,i])) - 1.17 
   endelse
endfor 
oplot, z, obyfe, color=3, linestyle=1

xyouts, 1.2, -0.7, 'M!Dhalo!N(z=0) = 10!E10.5!NM!D' + sun, charsize=1.5
multiplot

;-------------------------------------------------

halo_column = 170
restore, 'halo_template.sav'
o_data = read_ascii('set69/halos_o.out', template=stars_template)
fe_data =  read_ascii('set69/halos_Si.out', template=stars_template)
o = o_data.field001[halo_column,*]
fe = fe_data.field001[halo_column,*]
z = o_data.field001[0,*]
obyfe = fltarr(100)
for i = 0, 99 do begin 
   if (o[0,i] eq 0.0) then begin 
      obyfe[i] = -10.0
   endif else begin 
      obyfe[i] = alog10(abs(o[0,i]/fe[0,i])) - 1.17 
   endelse
endfor 
plot, z, obyfe, /xlog, xrange=[1,100], yrange=[-1,1], ytitle='[O/Si]'
oplot, z, obyfe, color=2

o_data = read_ascii('set67/halos_o.out', template=stars_template)
fe_data =  read_ascii('set67/halos_Si.out', template=stars_template)
o = o_data.field001[halo_column,*]
fe = fe_data.field001[halo_column,*]
z = o_data.field001[0,*]
obyfe = fltarr(100)
for i = 0, 99 do begin 
   if (o[0,i] eq 0.0) then begin 
      obyfe[i] = -10.0
   endif else begin 
      obyfe[i] = alog10(abs(o[0,i]/fe[0,i])) - 1.17 
   endelse
endfor 
oplot, z, obyfe, color=2, linestyle=3

o_data = read_ascii('set68/halos_o.out', template=stars_template)
fe_data =  read_ascii('set68/halos_Si.out', template=stars_template)
o = o_data.field001[halo_column,*]
fe = fe_data.field001[halo_column,*]
z = o_data.field001[0,*]
obyfe = fltarr(100)
for i = 0, 99 do begin 
   if (o[0,i] eq 0.0) then begin 
      obyfe[i] = -10.0
   endif else begin 
      obyfe[i] = alog10(abs(o[0,i]/fe[0,i])) - 1.17 
   endelse
endfor 
oplot, z, obyfe, color=2, linestyle=1

xyouts, 1.2, -0.7, 'M!Dhalo!N(z=0) = 10!E11.5!NM!D' + sun, charsize=1.5
multiplot, /reset 

;; device, /close_file
;; set_plot, 'X'

END
