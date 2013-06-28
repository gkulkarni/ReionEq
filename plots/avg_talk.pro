set_plot, 'ps'
device, filename='avg_obysi_talk.ps', xsize=7.0, ysize=3.0, $
        /inches, color=1, /HELVETICA
!P.font=0

;; window, xsize=500, ysize=300
;; Device, decomposed=0
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
!P.charsize = 1.5
!P.thick = 3.0
!P.charthick = 1

rs = fltarr(10)
cbyfe_avg = fltarr(10)
plot, rs, cbyfe_avg, yrange=[-1.0,1.0], ytitle='[O/Si]', $
      xrange=[2,6.5], xstyle=1, xtitle='redshift', /nodata,$
      xthick=3, ythick=3 

readcol, 'avg_data/cs2_51_obysi.dat', rs, cbyfe_avg, /silent 
oplot, rs, cbyfe_avg, thick=6
readcol, 'avg_data/cs2_49_obysi.dat', rs, cbyfe_avg, /silent 
oplot, rs, cbyfe_avg, color=2, thick=6
readcol, 'avg_data/cs2_50_obysi.dat', rs, cbyfe_avg, /silent 
oplot, rs, cbyfe_avg, color=3, thick=6

readcol, '../data/abr_cooke_sibyo.dat', z, s, sl, su, /silent
s = -s ; convert [Si/O] to [O/Si]
sl = -sl 
su = -su 
dy1 = su - s 
dy2 = s - sl 
plotsym, 0, 0.7, /FILL
oploterror, z, s, dy1, psym=8, /hibar
oploterror, z, s, dy2, psym=8, /lobar

readcol, '../data/abr_becker.dat', z, s, slerr, /silent
zt = z[2:*]
st = s[2:*]
st = -st ; convert [Si/O] to [O/Si]
oploterror, z, st, slerr, psym=8
zt = z[0:1]
st = s[0:1]
st = -st ; convert [Si/O] to [O/Si]
plotsym, 2, 2
oplot, zt, st, psym=8, thick=3

device, /close_file
set_plot, 'X'
