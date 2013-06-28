
; File: abevol.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; This procedure produces version of Figure 8.10 from
; Matteucci's old book in our model.  In other words, for a given halo
; mass, it plots the evolution of [O/Fe]. 

; lcolumn is column number in halos*.out files (starts from 0,
; corresponds to halo mass.) 

PRO abevol, lcolumn 

  window, xsize=1000, ysize=1000
  Device, decomposed=0
  TvLCT, 255, 0, 0, 2 
  TvLCT, 0, 127, 255, 3
  TvLCT, 255, 255, 0, 4
  !P.charsize = 2.0
  !P.thick = 1.0

  restore, 'halo_template.sav'
  fe_data = read_ascii('set33/halos_fe.out', template=stars_template)
  o_data = read_ascii('set33/halos_o.out', template=stars_template)

  z = fe_data.field001[0,*]
  fe = fe_data.field001[lcolumn,*]
  o = o_data.field001[lcolumn,*]

  obyfe = fltarr(100)

  for i = 0, 99 do begin 

     if (o[i] eq 0.0) then begin 
        obyfe[i] = 0.0
     endif else begin 
        obyfe[i] = alog10(abs(o[i]/fe[i])) - 0.91
     endelse

  endfor

  plot, z, obyfe, yrange=[-0.5,1.5], xrange=[50,0]

END


