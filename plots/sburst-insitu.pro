;; File: sburst-insitu.pro 
;;  Cre: 2014-01-09 

;; This is a simple IDL code to look at two particular files in
;; set183.  Probably not useful for anything else. 

window, xsize=750, ysize=750
!P.charsize = 2.0

openr, lun, 'set183/insitu', /get_lun
insitu_data = fltarr(271) 
readf, lun, insitu_data
insitu_data[0] = 0.0 
close, lun 
free_lun, lun 

openr, lun, 'set183/sburst', /get_lun
sburst_data = fltarr(271) 
readf, lun, sburst_data
sburst_data[0] = 0.0 
print, sburst_data
close, lun 
free_lun, lun 

sburst_data *= 0.1 
plot, insitu_data, /ylog, yrange=[1.0e-16,1.0e-12], xrange=[200,300]
oplot, sburst_data 

insitu_data += sburst_data 
oplot, insitu_data, thick=6
