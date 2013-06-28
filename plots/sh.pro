
; File: sh.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Plots mass and metallicity evolution of individual haloes from
; reion.  This plot was used in paper. 

;; set_plot, 'ps'
;; device, filename='sh.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

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
mmin_data = read_ascii('set71/mmin.out', template=mmintemplate) 
z = mmin_data.field1
mminc = mmin_data.field2
mminh = mmin_data.field3
mminav = mmin_data.field4

;; plot, z, mminav, /xlog, /ylog, xrange=[1,100], xtitle='!6z', ytitle='M!Dmin!N', linestyle=1
;; oplot, z, mminh, color=4, linestyle=5
;; oplot, z, mminc, color=3, linestyle=1

mminh=alog10(mminh)+10.0
mminav=alog10(mminav)+10.0
sol = "156B
sun = '!9' + string(sol) + '!X'
plot, z, mminh, /xlog, xrange=[1,100], linestyle=5, ytitle = 'log!D10!N(M/M!D' + sun + '!N)', yrange=[6,12]
oplot, z, mminav, linestyle=1

restore, 'halo_template.sav'
halos_data = read_ascii('set71/halos.out', template=stars_template)

halos = halos_data.field001[170,*]
halos=alog10(halos)+10.0
oplot, z, halos, color=2

halos = alog10(halos_data.field001[100,*])+10.0
oplot, z, halos, color=3

halos = alog10(halos_data.field001[73,*])+10.0
oplot, z, halos

xyouts, 30.0, 10.0, 'model 3', charsize=2.0, alignment=0.5

multiplot 

;---------------------------------------

;; zn_data = read_ascii('set71/halos_Zn.out', template=stars_template)
;; fe_data =  read_ascii('set71/halos_fe.out', template=stars_template)

;; zn = zn_data.field001[170,*]
;; fe = fe_data.field001[170,*]
;; znbyfe = fltarr(100)
;; for i = 0, 99 do begin 
;;    if (zn[0,i] eq 0.0) then begin 
;;       znbyfe[i] = -10.0
;;    endif else begin 
;;       znbyfe[i] = alog10(abs(zn[0,i]/fe[0,i])) + 3.07 
;;    endelse
;; endfor 
;; plot, z, znbyfe, /xlog, xrange=[1,100], yrange=[-1,0], xstyle=1
;; oplot, z, znbyfe, color=2

;; zn = zn_data.field001[100,*]
;; fe = fe_data.field001[100,*]
;; znbyfe = fltarr(100)
;; for i = 0, 99 do begin 
;;    if (zn[0,i] eq 0.0) then begin 
;;       znbyfe[i] = -10.0
;;    endif else begin 
;;       znbyfe[i] = alog10(abs(zn[0,i]/fe[0,i])) + 3.07 
;;    endelse
;; endfor 
;; oplot, z, znbyfe, color=3

;; zn = zn_data.field001[100,*]
;; fe = fe_data.field001[100,*]
;; znbyfe = fltarr(100)
;; for i = 0, 99 do begin 
;;    if (zn[0,i] eq 0.0) then begin 
;;       znbyfe[i] = -10.0
;;    endif else begin 
;;       znbyfe[i] = alog10(abs(zn[0,i]/fe[0,i])) + 3.07 
;;    endelse
;; endfor 
;; oplot, z, znbyfe, color=4

;---------------------------------------

;; solar = -1.87 ; O/H 
;; o_data = read_ascii('set22/halos_o.out', template=stars_template)
;; gas_data =  read_ascii('set22/halos_gas.out', template=stars_template)

;; o = o_data.field001[200,*]
;; gas = gas_data.field001[200,*]
;; obyh = fltarr(100)
;; for i = 0, 99 do begin 
;;    mH = gas[0,i]*0.71 
;;    if (o[0,i] eq 0.0) then begin 
;;       obyh[i] = -10.0
;;    endif else begin 
;;       obyh[i] = alog10(abs(o[0,i]/mH)) - solar 
;;    endelse
;; endfor 
;; plot, z, obyh, /xlog, xrange=[1,100], yrange=[-8,0], ytitle='[O/H]'
;; oplot, z, obyh, color=2

;; o = o_data.field001[150,*]
;; gas = gas_data.field001[150,*]
;; for i = 0, 99 do begin 
;;    mH = gas[0,i]*0.71 
;;    if (o[0,i] eq 0.0) then begin 
;;       obyh[i] = -10.0
;;    endif else begin 
;;       obyh[i] = alog10(abs(o[0,i]/mH)) - solar 
;;    endelse
;; endfor 
;; oplot, z, obyh, color=3

;; o = o_data.field001[100,*]
;; gas = gas_data.field001[100,*]
;; for i = 0, 99 do begin 
;;    mH = gas[0,i]*0.71 
;;    if (o[0,i] eq 0.0) then begin 
;;       obyh[i] = -10.0
;;    endif else begin 
;;       obyh[i] = alog10(abs(o[0,i]/mH)) - solar 
;;    endelse
;; endfor 
;; oplot, z, obyh

;----------------------------------------------------------------

;; solar = -2.78 ; O/H 
;; o_data = read_ascii('set22/halos_fe.out', template=stars_template)
;; gas_data =  read_ascii('set22/halos_gas.out', template=stars_template)

;; o = o_data.field001[170,*]
;; gas = gas_data.field001[170,*]
;; obyh = fltarr(100)
;; for i = 0, 99 do begin 
;;    mH = gas[0,i]*0.71 
;;    if (o[0,i] eq 0.0) then begin 
;;       obyh[i] = -10.0
;;    endif else begin 
;;       obyh[i] = alog10(abs(o[0,i]/mH)) - solar 
;;    endelse
;; endfor 
;; plot, z, obyh, /xlog, xrange=[1,100], yrange=[-4,-1.5], ytitle='[Fe/H]', xstyle=1
;; oplot, z, obyh, color=2

;; o = o_data.field001[100,*]
;; gas = gas_data.field001[100,*]
;; for i = 0, 99 do begin 
;;    mH = gas[0,i]*0.71 
;;    if (o[0,i] eq 0.0) then begin 
;;       obyh[i] = -10.0
;;    endif else begin 
;;       obyh[i] = alog10(abs(o[0,i]/mH)) - solar 
;;    endelse
;; endfor 
;; oplot, z, obyh, color=3

;; o = o_data.field001[100,*]
;; gas = gas_data.field001[100,*]
;; for i = 0, 99 do begin 
;;    mH = gas[0,i]*0.71 
;;    if (o[0,i] eq 0.0) then begin 
;;       obyh[i] = -10.0
;;    endif else begin 
;;       obyh[i] = alog10(abs(o[0,i]/mH)) - solar 
;;    endelse
;; endfor 
;; oplot, z, obyh

;----------------------------------------------------------------

solar = 1.17
o_data = read_ascii('set71/halos_o.out', template=stars_template)
si_data =  read_ascii('set71/halos_Si.out', template=stars_template)

o = o_data.field001[73,*]
print, o 
si = si_data.field001[73,*]
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
; oplot, z, obysi, color=2


o = o_data.field001[100,*]
si = si_data.field001[100,*]
obysi = fltarr(len)
for i = 0, len-1 do begin 
   if (o[0,i] eq 0.0) then begin 
      obysi[i] = -10.0
   endif else begin 
      obysi[i] = alog10(abs(o[0,i]/si[0,i])) - solar 
   endelse
endfor 
oplot, z, obysi, color=3

o = o_data.field001[170,*]
si = si_data.field001[170,*]
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

END
