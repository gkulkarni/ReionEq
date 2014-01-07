
; File: sh.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Plots mass and metallicity evolution of individual haloes from
; reion.  This plot was used in paper. 

;; set_plot, 'ps'
;; device, filename='sh.ps', xsize=7.0, ysize=7.0, /inches, color=1,
;; yoffset=1.0

PRO sh, opt 

  window, xsize=700, ysize=700
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
  mmin_data = read_ascii('set159/mmin.out', template=mmintemplate) 
  z = mmin_data.field1
  mminc = mmin_data.field2
  mminh = mmin_data.field3
  mminav = mmin_data.field4

  mminh=alog10(mminh)+10.0
  mminav=alog10(mminav)+10.0
  plot, z, mminh, /xlog, xrange=[1,100], linestyle=5, ytitle = 'log!D10!N(M/M!D!9n!X!N)', yrange=[6,12]
  oplot, z, mminav, linestyle=1

  restore, 'halo_template.sav'
  halos_data = read_ascii('set159/halos.out', template=stars_template)

  halos = halos_data.field001[250,*]
  halos=alog10(halos)+10.0
  oplot, z, halos, color=2

  halos = alog10(halos_data.field001[140,*])+10.0
  oplot, z, halos, color=3

  halos = alog10(halos_data.field001[100,*])+10.0
  oplot, z, halos

  xyouts, 30.0, 10.0, 'model 3', charsize=2.0, alignment=0.5

  multiplot 

  case opt of 

     1: begin 
        solar = 1.17
        o_data = read_ascii('set126/halos_o.out', template=stars_template)
        si_data =  read_ascii('set126/halos_Si.out', template=stars_template)

        o = o_data.field001[100,*]
        si = si_data.field001[100,*]
        len = size(o, /n_elements)
        obysi = fltarr(len)
        for i = 0, len-1 do begin 
           if (o[0,i] eq 0.0) then begin 
              obysi[i] = -10.0
           endif else begin 
              obysi[i] = alog10(abs(o[0,i]/si[0,i])) - solar 
           endelse
        endfor 
        plot, z, obysi, /xlog, xrange=[1,100], yrange=[-1,1], ytitle='[O/Si]', xstyle=1

        o = o_data.field001[140,*]
        si = si_data.field001[140,*]
        obysi = fltarr(len)
        for i = 0, len-1 do begin 
           if (o[0,i] eq 0.0) then begin 
              obysi[i] = -10.0
           endif else begin 
              obysi[i] = alog10(abs(o[0,i]/si[0,i])) - solar 
           endelse
        endfor 
        oplot, z, obysi, color=3

        o = o_data.field001[250,*]
        si = si_data.field001[250,*]
        obysi = fltarr(len)
        for i = 0, len-1 do begin 
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
     end

     2: begin 
        solar = 0.41
        c_data = read_ascii('set125/halos_c.out', template=stars_template)
        fe_data =  read_ascii('set125/halos_fe.out', template=stars_template)

        c = c_data.field001[100,*]
        fe = fe_data.field001[100,*]
        len = size(c, /n_elements)
        cbyfe = fltarr(len)
        for i = 0, len-1 do begin 
           if (c[0,i] eq 0.0) then begin 
              cbyfe[i] = -10.0
           endif else begin 
              cbyfe[i] = alog10(abs(c[0,i]/fe[0,i])) - solar 
           endelse
        endfor 
        plot, z, cbyfe, /xlog, xrange=[1,100], yrange=[-1,1], ytitle='[C/Fe]', xstyle=1

        c = c_data.field001[140,*]
        fe = fe_data.field001[140,*]
        cbyfe = fltarr(len)
        for i = 0, len-1 do begin 
           if (c[0,i] eq 0.0) then begin 
              cbyfe[i] = -10.0
           endif else begin 
              cbyfe[i] = alog10(abs(c[0,i]/fe[0,i])) - solar 
           endelse
        endfor 
        oplot, z, cbyfe, color=3

        c = c_data.field001[250,*]
        fe = fe_data.field001[250,*]
        cbyfe = fltarr(len)
        for i = 0, len-1 do begin 
           if (c[0,i] eq 0.0) then begin 
              cbyfe[i] = -10.0
           endif else begin 
              cbyfe[i] = alog10(abs(c[0,i]/fe[0,i])) - solar 
           endelse
        endfor 
        oplot, z, cbyfe, color=2
        vline, 6.0, linestyle=2

        multiplot, /reset 
        close, /all 
     end
  endcase

END
