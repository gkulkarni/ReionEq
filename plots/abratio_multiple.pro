
PRO abratio, row, redshift, set1, set2, set3

window, xsize=2500, ysize=1000
!P.multi = [0, 5, 2, 0, 1]
!P.charsize = 3.0
!y.omargin = [2,6]
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3

xsize = 16.0
margin = 0.12
wall = 0.12

a = xsize/32.0 - (margin + wall) 
b = a 
ysize = (margin + b + wall + b + wall)*16.0

; row = 95
redshift_label = 'z=2'
redshift_label = 'z='+strtrim(uint(redshift),2)

file_prefix1 = 'set'+strtrim(uint(set1),2)
file_prefix2 = 'set'+strtrim(uint(set2),2)
file_prefix3 = 'set'+strtrim(uint(set3),2)

lcolumn = 95
restore, 'halo_template.sav'

c_data1 = read_ascii(file_prefix1+'/halos_c.out', template=stars_template)
fe_data1 = read_ascii(file_prefix1+'/halos_fe.out', template=stars_template)
o_data1 = read_ascii(file_prefix1+'/halos_o.out', template=stars_template)
n_data1 = read_ascii(file_prefix1+'/halos_N.out', template=stars_template)
si_data1 =  read_ascii(file_prefix1+'/halos_Si.out', template=stars_template)
zn_data1 = read_ascii(file_prefix1+'/halos_Zn.out', template=stars_template)
mg_data1 = read_ascii(file_prefix1+'/halos_Mg.out', template=stars_template)
gas_data1 =  read_ascii(file_prefix1+'/halos_gas.out', template=stars_template)
halos_data1 = read_ascii(file_prefix1+'/halos.out', template=stars_template)

c1 = c_data1.field001[lcolumn:*,row]
fe1 = fe_data1.field001[lcolumn:*,row]
o1 = o_data1.field001[lcolumn:*,row]
n1 = n_data1.field001[lcolumn:*,row]
si1 = si_data1.field001[lcolumn:*,row]
zn1 = zn_data1.field001[lcolumn:*,row]
mg1 = mg_data1.field001[lcolumn:*,row]
gas1 = gas_data1.field001[lcolumn:*,row]
halos1 = halos_data1.field001[lcolumn:*,row]

c_data2 = read_ascii(file_prefix2+'/halos_c.out', template=stars_template)
fe_data2 = read_ascii(file_prefix2+'/halos_fe.out', template=stars_template)
o_data2 = read_ascii(file_prefix2+'/halos_o.out', template=stars_template)
n_data2 = read_ascii(file_prefix2+'/halos_N.out', template=stars_template)
si_data2 = read_ascii(file_prefix2+'/halos_Si.out', template=stars_template)
zn_data2 = read_ascii(file_prefix2+'/halos_Zn.out', template=stars_template)
mg_data2 = read_ascii(file_prefix2+'/halos_Mg.out', template=stars_template)
gas_data2 = read_ascii(file_prefix2+'/halos_gas.out', template=stars_template)
halos_data2 = read_ascii(file_prefix2+'/halos.out', template=stars_template)

c2 = c_data2.field001[lcolumn:*,row]
fe2 = fe_data2.field001[lcolumn:*,row]
o2 = o_data2.field001[lcolumn:*,row]
n2 = n_data2.field001[lcolumn:*,row]
si2 = si_data2.field001[lcolumn:*,row]
zn2 = zn_data2.field001[lcolumn:*,row]
mg2 = mg_data2.field001[lcolumn:*,row]
gas2 = gas_data2.field001[lcolumn:*,row]
halos2 = halos_data2.field001[lcolumn:*,row]

c_data3 = read_ascii(file_prefix3+'/halos_c.out', template=stars_template)
fe_data3 = read_ascii(file_prefix3+'/halos_fe.out', template=stars_template)
o_data3 = read_ascii(file_prefix3+'/halos_o.out', template=stars_template)
n_data3 = read_ascii(file_prefix3+'/halos_N.out', template=stars_template)
si_data3 = read_ascii(file_prefix3+'/halos_Si.out', template=stars_template)
zn_data3 = read_ascii(file_prefix3+'/halos_Zn.out', template=stars_template)
mg_data3 = read_ascii(file_prefix3+'/halos_Mg.out', template=stars_template)
gas_data3 = read_ascii(file_prefix3+'/halos_gas.out', template=stars_template)
halos_data3 = read_ascii(file_prefix3+'/halos.out', template=stars_template)
halos3 = halos_data3.field001[lcolumn:*,row]

c3 = c_data3.field001[lcolumn:*,row]
fe3 = fe_data3.field001[lcolumn:*,row]
o3 = o_data3.field001[lcolumn:*,row]
n3 = n_data3.field001[lcolumn:*,row]
si3 = si_data3.field001[lcolumn:*,row]
zn3 = zn_data3.field001[lcolumn:*,row]
mg3 = mg_data3.field001[lcolumn:*,row]
gas3 = gas_data3.field001[lcolumn:*,row]

array_size = size(c3)
array_length = array_size[1] 

cbyfe1 = fltarr(array_length)
cbyfe2 = cbyfe1
cbyfe3 = cbyfe1

febyh1 = cbyfe1
febyh2 = cbyfe1
febyh3 = cbyfe1

znbyfe1 = cbyfe1
znbyfe2 = cbyfe1
znbyfe3 = cbyfe1

obyfe1 = cbyfe1
obyfe2 = cbyfe1
obyfe3 = cbyfe1

sibyfe1 = cbyfe1
sibyfe2 = cbyfe1
sibyfe3 = cbyfe1

obyh1 = cbyfe1
obyh2 = cbyfe1
obyh3 = cbyfe1

cbyo1 = cbyfe1
cbyo2 = cbyfe1
cbyo3 = cbyfe1

nbysi1 = cbyfe1
nbysi2 = cbyfe1
nbysi3 = cbyfe1

sibyh1 = cbyfe1
sibyh2 = cbyfe1
sibyh3 = cbyfe1

obysi1 = cbyfe1
obysi2 = cbyfe1
obysi3 = cbyfe1

mgbyfe1 = cbyfe1 
mgbyfe2 = cbyfe1 
mgbyfe3 = cbyfe1 

nbyo1 = cbyfe1 
nbyo2 = cbyfe1 
nbyo3 = cbyfe1 


for i = 0, array_length-1 do begin 
   
   mH1 = gas1[i]*0.71
   
   if (fe1[i] eq 0.0) then begin 
      febyh1[i] = 0.0
   endif else begin 
      febyh1[i] = alog10(abs(fe1[i]/mH1)) + 2.78 
   endelse

   if (c1[i] eq 0.0) then begin 
      cbyfe1[i] = 0.0
   endif else begin 
      cbyfe1[i] = alog10(abs(c1[i]/fe1[i])) - 0.41
   endelse

   if (zn1[i] eq 0.0) then begin 
      znbyfe1[i] = 0.0
   endif else begin 
      znbyfe1[i] = alog10(abs(zn1[i]/fe1[i])) + 3.07
   endelse

   if (n1[i] eq 0.0) then begin 
      nbysi1[i] = 0.0
   endif else begin 
      nbysi1[i] = alog10(abs(n1[i]/si1[i])) - 0.22
   endelse

   if (si1[i] eq 0.0) then begin 
      sibyh1[i] = 0.0
   endif else begin 
      sibyh1[i] = alog10(abs(si1[i]/mH1)) + 3.03
   endelse

   if (c1[i] eq 0.0) then begin 
      cbyo1[i] = 0.0
   endif else begin 
      cbyo1[i] = alog10(abs(c1[i]/o1[i])) + 0.5
   endelse

   if (o1[i] eq 0.0) then begin 
      obyh1[i] = 0.0
   endif else begin 
      obyh1[i] = alog10(abs(o1[i]/mH1)) + 1.87
   endelse

   if (o1[i] eq 0.0) then begin 
      obyfe1[i] = 0.0
   endif else begin 
      obyfe1[i] = alog10(abs(o1[i]/fe1[i])) - 0.91
   endelse

   if (si1[i] eq 0.0) then begin 
      sibyfe1[i] = 0.0
   endif else begin 
      sibyfe1[i] = alog10(abs(si1[i]/fe1[i])) + 0.25 
   endelse

   if (si1[i] eq 0.0) then begin 
      obysi1[i] = 0.0
   endif else begin 
      obysi1[i] = alog10(abs(o1[i]/si1[i])) - 1.17 
   endelse

   if (mg1[i] eq 0.0) then begin 
      mgbyfe1[i] = 0.0
   endif else begin 
      mgbyfe1[i] = alog10(abs(mg1[i]/fe1[i])) - 0.44
   endelse

   if (n1[i] eq 0.0) then begin 
      nbyo1[i] = 0.0
   endif else begin 
      nbyo1[i] = alog10(abs(n1[i]/o1[i])) + 0.94 
   endelse

endfor 

for i = 0, array_length-1 do begin 
   
   mH2 = gas2[i]*0.71
   
   if (fe2[i] eq 0.0) then begin 
      febyh2[i] = 0.0
   endif else begin 
      febyh2[i] = alog10(abs(fe2[i]/mH2)) + 2.78 
   endelse

   if (c2[i] eq 0.0) then begin 
      cbyfe2[i] = 0.0
   endif else begin 
      cbyfe2[i] = alog10(abs(c2[i]/fe2[i])) - 0.41
   endelse

   if (zn2[i] eq 0.0) then begin 
      znbyfe2[i] = 0.0
   endif else begin 
      znbyfe2[i] = alog10(abs(zn2[i]/fe2[i])) + 3.07
   endelse

   if (n2[i] eq 0.0) then begin 
      nbysi2[i] = 0.0
   endif else begin 
      nbysi2[i] = alog10(abs(n2[i]/si2[i])) - 0.22
   endelse

   if (si2[i] eq 0.0) then begin 
      sibyh2[i] = 0.0
   endif else begin 
      sibyh2[i] = alog10(abs(si2[i]/mH2)) + 3.03
   endelse

   if (c2[i] eq 0.0) then begin 
      cbyo2[i] = 0.0
   endif else begin 
      cbyo2[i] = alog10(abs(c2[i]/o2[i])) + 0.5
   endelse

   if (o2[i] eq 0.0) then begin 
      obyh2[i] = 0.0
   endif else begin 
      obyh2[i] = alog10(abs(o2[i]/mH2)) + 1.87
   endelse

   if (o2[i] eq 0.0) then begin 
      obyfe2[i] = 0.0
   endif else begin 
      obyfe2[i] = alog10(abs(o2[i]/fe2[i])) - 0.91
   endelse

   if (si2[i] eq 0.0) then begin 
      sibyfe2[i] = 0.0
   endif else begin 
      sibyfe2[i] = alog10(abs(si2[i]/fe2[i])) + 0.25 
   endelse

   if (si2[i] eq 0.0) then begin 
      obysi2[i] = 0.0
   endif else begin 
      obysi2[i] = alog10(abs(o2[i]/si2[i])) - 1.17
   endelse

   if (mg2[i] eq 0.0) then begin 
      mgbyfe2[i] = 0.0
   endif else begin 
      mgbyfe2[i] = alog10(abs(mg2[i]/fe2[i])) - 0.44
   endelse

   if (n2[i] eq 0.0) then begin 
      nbyo2[i] = 0.0
   endif else begin 
      nbyo2[i] = alog10(abs(n2[i]/o2[i])) + 0.94 
   endelse

endfor 

for i = 0, array_length-1 do begin 
   
   mH3 = gas3[i]*0.71
   
   if (fe3[i] eq 0.0) then begin 
      febyh3[i] = 0.0
   endif else begin 
      febyh3[i] = alog10(abs(fe3[i]/mH3)) + 2.78 
   endelse

   if (c3[i] eq 0.0) then begin 
      cbyfe3[i] = 0.0
   endif else begin 
      cbyfe3[i] = alog10(abs(c3[i]/fe3[i])) - 0.41
   endelse

   if (zn3[i] eq 0.0) then begin 
      znbyfe3[i] = 0.0
   endif else begin 
      znbyfe3[i] = alog10(abs(zn3[i]/fe3[i])) + 3.07
   endelse

   if (n3[i] eq 0.0) then begin 
      nbysi3[i] = 0.0
   endif else begin 
      nbysi3[i] = alog10(abs(n3[i]/si3[i])) - 0.22
   endelse

   if (si3[i] eq 0.0) then begin 
      sibyh3[i] = 0.0
   endif else begin 
      sibyh3[i] = alog10(abs(si3[i]/mH3)) + 3.03
   endelse

   if (c3[i] eq 0.0) then begin 
      cbyo3[i] = 0.0
   endif else begin 
      cbyo3[i] = alog10(abs(c3[i]/o3[i])) + 0.5
   endelse

   if (o3[i] eq 0.0) then begin 
      obyh3[i] = 0.0
   endif else begin 
      obyh3[i] = alog10(abs(o3[i]/mH3)) + 1.87
   endelse

   if (o3[i] eq 0.0) then begin 
      obyfe3[i] = 0.0
   endif else begin 
      obyfe3[i] = alog10(abs(o3[i]/fe3[i])) - 0.91
   endelse

   if (si3[i] eq 0.0) then begin 
      sibyfe3[i] = 0.0
   endif else begin 
      sibyfe3[i] = alog10(abs(si3[i]/fe3[i])) + 0.25 
   endelse

   if (si3[i] eq 0.0) then begin 
      obysi3[i] = 0.0
   endif else begin 
      obysi3[i] = alog10(abs(o3[i]/si3[i])) - 1.17
   endelse

   if (mg3[i] eq 0.0) then begin 
      mgbyfe3[i] = 0.0
   endif else begin 
      mgbyfe3[i] = alog10(abs(mg3[i]/fe3[i])) - 0.44
   endelse

   if (n3[i] eq 0.0) then begin 
      nbyo3[i] = 0.0
   endif else begin 
      nbyo3[i] = alog10(abs(n3[i]/o3[i])) + 0.94 
   endelse

endfor 

f2 = smooth(febyh1, 30, /edge_truncate)
f3 = smooth(f2, 30, /edge_truncate)
c2 = smooth(cbyfe1, 80, /edge_truncate)
c3 = smooth(c2, 20, /edge_truncate)
plot, f3, c3, xtitle='!6[Fe/H]', ytitle='[C/Fe]', yrange=[-1,1], xrange=[-2.5,-1]

f3 = smooth(febyh2, 30, /edge_truncate)
c2 = smooth(cbyfe2, 80, /edge_truncate)
c3 = smooth(c2, 20, /edge_truncate)
oplot, f3, c3, color=2

f3 = smooth(febyh3, 30, /edge_truncate)
c2 = smooth(cbyfe3, 80, /edge_truncate)
c3 = smooth(c2, 20, /edge_truncate)
oplot, f3, c3, color=3

;------------------------------------------------

f2 = smooth(febyh1, 30, /edge_truncate)
f3 = smooth(f2, 30, /edge_truncate)
zn2 = smooth(znbyfe1, 30, /edge_truncate)
zn3=smooth(zn2, 30, /edge_truncate)
plot, f3, zn3, xtitle='[Fe/H]', ytitle='[Zn/Fe]', yrange=[-1,1], xrange=[-2.5,-1]

f3 = smooth(febyh2, 30, /edge_truncate)
zn2 = smooth(znbyfe2, 30, /edge_truncate)
zn3=smooth(zn2, 30, /edge_truncate)
oplot, f3, zn3, color=2

f3 = smooth(febyh3, 30, /edge_truncate)
zn2 = smooth(znbyfe3, 30, /edge_truncate)
zn3=smooth(zn2, 30, /edge_truncate)
oplot, f3, zn3, color=3

;------------------------------------------------

f3 = smooth(febyh1, 30, /edge_truncate)
mg2 = smooth(mgbyfe1, 30, /edge_truncate)
mg3=smooth(mg2, 30, /edge_truncate)
plot, f3, mg3, ytitle='[Mg/Fe]', xtitle='[Fe/H]', yrange=[-1,1], xrange=[-2.5,-1]

f3 = smooth(febyh2, 30, /edge_truncate)
mg2 = smooth(mgbyfe2, 30, /edge_truncate)
mg3=smooth(mg2, 30, /edge_truncate)
oplot, f3, mg3, color=2

f3 = smooth(febyh3, 30, /edge_truncate)
mg2 = smooth(mgbyfe3, 30, /edge_truncate)
mg3=smooth(mg2, 30, /edge_truncate)
oplot, f3, mg3, color=3

;------------------------------------------------

f3 = smooth(febyh1, 30, /edge_truncate)
si2 = smooth(sibyfe1, 30, /edge_truncate)
si3=smooth(si2, 30, /edge_truncate)
plot, f3, si3, ytitle='[Si/Fe]', xtitle='[Fe/H]', yrange=[-1,1], xrange=[-2.5,-1]

f3 = smooth(febyh2, 30, /edge_truncate)
si2 = smooth(sibyfe2, 30, /edge_truncate)
si3=smooth(si2, 30, /edge_truncate)
oplot, f3, si3, color=2

f3 = smooth(febyh3, 30, /edge_truncate)
si2 = smooth(sibyfe3, 30, /edge_truncate)
si3=smooth(si2, 30, /edge_truncate)
oplot, f3, si3, color=3

;------------------------------------------------

f2 = smooth(febyh1, 30, /edge_truncate)
f3 = smooth(f2, 30, /edge_truncate)
o2 = smooth(obyfe1, 30, /edge_truncate)
o3=smooth(o2, 30, /edge_truncate)
plot, f3, o3, xtitle='[Fe/H]', ytitle='[O/Fe]', yrange=[-1,1], xrange=[-2.5,-1]

;; opick = fltarr(1)
;; fpick = fltarr(1)
;; opick[0] = o3[hpick]
;; fpick[0] = f3[hpick]
;; plotsym, 0, 1, /FILL
;; oplot, fpick, opick, psym=8

f3 = smooth(febyh2, 30, /edge_truncate)
o2 = smooth(obyfe2, 30, /edge_truncate)
o3=smooth(o2, 30, /edge_truncate)
oplot, f3, o3, color=2

;; opick[0] = o3[hpick]
;; fpick[0] = f3[hpick]
;; plotsym, 0, 1, /FILL
;; oplot, fpick, opick, psym=8, color=2

f3 = smooth(febyh3, 30, /edge_truncate)
o2 = smooth(obyfe3, 30, /edge_truncate)
o3=smooth(o2, 30, /edge_truncate)
oplot, f3, o3, color=3

;; opick[0] = o3[hpick]
;; fpick[0] = f3[hpick]
;; plotsym, 0, 1, /FILL
;; oplot, fpick, opick, psym=8, color=3

;------------------------------------------------

o2 = smooth(obyh1, 30, /edge_truncate)
o3 = smooth(o2, 30, /edge_truncate) 
c2 = smooth(obysi1, 30, /edge_truncate)
c3=smooth(c2, 30, /edge_truncate)
plot, o3, c3, xtitle='[O/H]', ytitle='[O/Si]', yrange=[-1,1], xrange=[-1.8,-0.8]

o2 = smooth(obyh2, 30, /edge_truncate)
o3 = smooth(o2, 30, /edge_truncate) 
c2 = smooth(obysi2, 30, /edge_truncate)
c3=smooth(c2, 30, /edge_truncate)
oplot, o3, c3, color=2

o2 = smooth(obyh3, 30, /edge_truncate)
c2 = smooth(obysi3, 30, /edge_truncate)
c3=smooth(c2, 30, /edge_truncate)
oplot, o2, c3, color=3

;------------------------------------------------

o2 = smooth(obyh1, 30, /edge_truncate)
o3 = smooth(o2, 30, /edge_truncate) 
c2 = smooth(cbyo1, 30, /edge_truncate)
c3=smooth(c2, 30, /edge_truncate)
plot, o3, c3, xtitle='[O/H]', ytitle='[C/O]', yrange=[-1,1], xrange=[-1.8,-0.8]

o2 = smooth(obyh2, 30, /edge_truncate)
o3 = smooth(o2, 30, /edge_truncate) 
c2 = smooth(cbyo2, 30, /edge_truncate)
c3=smooth(c2, 30, /edge_truncate)
oplot, o3, c3, color=2

o2 = smooth(obyh3, 30, /edge_truncate)
c2 = smooth(cbyo3, 30, /edge_truncate)
c3=smooth(c2, 30, /edge_truncate)
oplot, o2, c3, color=3

;------------------------------------------------

si2 = smooth(sibyh1, 30, /edge_truncate)
n2 = smooth(nbysi1, 30, /edge_truncate)
plot, si2, n2, xtitle='[Si/H]', ytitle='[N/Si]', yrange=[-3,-1], xrange=[-1.8,-1.0]

si2 = smooth(sibyh2, 30, /edge_truncate)
n2 = smooth(nbysi2, 30, /edge_truncate)
oplot, si2, n2, color=2

si2 = smooth(sibyh3, 30, /edge_truncate)
n2 = smooth(nbysi3, 30, /edge_truncate)
oplot, si2, n2, color=3


;------------------------------------------------

n2 = smooth(nbyo1, 30, /edge_truncate)
o2 = smooth(obyh1, 30, /edge_truncate)
plot, o2, n2, xtitle='[O/H]', ytitle='[N/O]', yrange=[-3.5,-2.5], xrange=[-1.8,-0.8]

n2 = smooth(nbyo2, 30, /edge_truncate)
o2 = smooth(obyh2, 30, /edge_truncate)
oplot, o2, n2, color=2

n2 = smooth(nbyo3, 30, /edge_truncate)
o2 = smooth(obyh3, 30, /edge_truncate)
oplot, o2, n2, color=3

;------------------------------------------------

XYOuts, 0.5, 0.95, ALIGNMENT=0.5, CHARSIZE=4.0, /NORMAL, redshift_label


END
