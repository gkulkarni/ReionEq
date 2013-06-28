
; File: pop3.pro 
;  Cre: 2012
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Plots redshift evolution of number of Pop III mass bins in
; reion. Compares it with redshift evolution of mmin.f90.

restore, 'halo_template.sav'
pop_data = read_ascii('set69/strpop.out', template=stars_template)
p3_count = intarr(100)
for row = 0, 99 do begin
   pop = pop_data.field001[1:*, row] 
   pop = pop - 2 
   index = where(pop, count) 
;   print, row, count 
   p3_count[row] = count 
endfor 

restore, 'mmintemplate.sav'
mmin_data = read_ascii('set69/mmin.out', template=mmintemplate) 
z = mmin_data.field1

plot, z, p3_count, psym=-6, /xlog, xrange=[1,100], thick=2

pop_data = read_ascii('set67/strpop.out', template=stars_template)
p3_count = intarr(100)
for row = 0, 99 do begin
   pop = pop_data.field001[1:*, row] 
   pop = pop - 2 
   index = where(pop, count) 
;   print, row, count 
   p3_count[row] = count 
endfor 
mmin_data = read_ascii('set67/mmin.out', template=mmintemplate) 
z = mmin_data.field1

oplot, z, p3_count, color=2, psym=-5

pop_data = read_ascii('set68/strpop.out', template=stars_template)
p3_count = intarr(100)
for row = 0, 99 do begin
   pop = pop_data.field001[1:*, row] 
   pop = pop - 2 
   index = where(pop, count) 
;   print, row, count 
   p3_count[row] = count 
endfor 

mmin_data = read_ascii('set68/mmin.out', template=mmintemplate) 
z = mmin_data.field1

oplot, z, p3_count, color=3, psym=-6

END
