;; set_plot, 'ps'
;; device, filename='zism.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0
window, xsize=1000, ysize=1000
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 177, 218, 226, 5 
!P.charsize = 2.0
!P.thick = 3.0
!P.charthick = 3 

restore, 'abundancetemplate.sav'
abundance_data = read_ascii('set49/abundances.out', template=abundancetemplate)
; abundance_data = read_ascii('set5/abundances.out', template=abundancetemplate)
z=abundance_data.field1 
zism=abundance_data.field2 
sol = "156B
sun = '!9' + string(sol) + '!X'
plot, z, zism, /ylog, /xlog, xrange=[1,50], xstyle=1, yrange=[1.0e-10,1.0e-1], xtitle='!6z', ytitle='Z/Z!D' + sun

openr, lun, '../data/prochaska03.dat', /get_lun
obsdata_zism = fltarr(2, 109)
readf, lun, obsdata_zism
z_obs = obsdata_zism[0,*]
zism_obs = obsdata_zism[1,*]
oplot, z_obs, zism_obs, psym=7, color=2

restore, 'fractemplate.sav'
fracs_data = read_ascii('set40/fracs.out', template=fractemplate) 
z = fracs_data.field1
xigm = fracs_data.field5
zigm = xigm / 0.02
xism = fracs_data.field9
zism2 = xism / 0.02 
;oplot, z, zigm, linestyle=5

;; abundance_data = read_ascii('set15/abundances.out', template=abundancetemplate)
;; z=abundance_data.field1 
;; zism=abundance_data.field2 
;; ;oplot, z, zism, color=3

;; fracs_data = read_ascii('set15/fracs.out', template=fractemplate) 
;; z = fracs_data.field1
;; xigm = fracs_data.field5
;; zigm = xigm / 0.02
;; ;oplot, z, zigm, linestyle=5, color=3

;; abundance_data = read_ascii('set34/abundances.out', template=abundancetemplate)
;; z=abundance_data.field1
;; zism=abundance_data.field2
;; oplot, z, zism, color=2

;; abundance_data = read_ascii('set35/abundances.out', template=abundancetemplate)
;; z=abundance_data.field1
;; zism=abundance_data.field2
;; oplot, z, zism, color=3

;; device, /close_file
;; set_plot, 'X'
