
; File: avg_multi.pro 
;  Cre: 2012 
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; This code plots the observations of Becker et al. 2011.  It then
; overplots our theoretical calculation of the same quantities.  This
; plot was used in the paper.

;; set_plot, 'ps'
;; device, filename='avg_multi2.ps', xsize=7.0, ysize=5.0, /inches, color=1, yoffset=1.0

window, xsize=1000, ysize=600
Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
!P.charsize = 1.5
!P.thick = 1.0
!P.charthick = 1

erase
multiplot, [1,2], mXtitle='!6z', mXtitsize='1.5', mytitsize='1.5', $
           mYtitoffset=1.5, xtickformat='(I1)', mXtitoffset=1.5, /doyaxis
@avg
plotsym, 0, 0.7, /FILL
readcol, '../data/abr_cooke_cbyfe.dat', z, c, cl, cu 
oplot, z, c, psym=8 
dy1 = cu - c 
dy2 = c - cl 
oploterror, z, c, dy1, psym=8, /hibar
oploterror, z, c, dy2, psym=8, /lobar
readcol, '../data/abr_becker_cbyfe.dat', z, c, cerr 
zt = fltarr(4) 
ct = fltarr(4) 
zt2 = zt 
ct2 = ct 
cerrt = fltarr(4) 
j = 0 
k = 0 
for i = 0, 7 do begin 
   if cerr[i] ne 0.0 then begin 
      zt[j] = z[i] 
      ct[j] = c[i] 
      cerrt[j] = cerr[i] 
      j = j+1 
   endif else begin 
      zt2[k] = z[i] 
      ct2[k] = c[i] 
      k = k+1 
   endelse
endfor
oploterror, zt, ct, cerrt, psym=8
plotsym, 2, 2
oplot, zt2, ct2, psym=8 
multiplot
@avg2
readcol, '../data/abr_cooke_sibyo.dat', z, s, sl, su 
s = -s ; convert [Si/O] to [O/Si]
sl = -sl 
su = -su 
dy1 = su - s 
dy2 = s - sl 
plotsym, 0, 0.7, /FILL
oploterror, z, s, dy1, psym=8, /hibar
oploterror, z, s, dy2, psym=8, /lobar
readcol, '../data/abr_becker.dat', z, s, slerr 
zt = z[2:*]
st = s[2:*]
st = -st ; convert [Si/O] to [O/Si]
oploterror, z, st, slerr, psym=8
zt = z[0:1]
st = s[0:1]
st = -st ; convert [Si/O] to [O/Si]
plotsym, 2, 2
oplot, zt, st, psym=8
;multiplot
;@avg3
multiplot, /reset 
close, /all

;; device, /close_file
;; set_plot, 'X'

END

