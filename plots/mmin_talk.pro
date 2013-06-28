; File: mmin_talk.pro 
;  Cre: 2013-02-06
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $)

set_plot, 'ps'
device, filename='mmin_talk.ps', xsize=7.0, ysize=7.0, $
        /inches, color=1, /HELVETICA, yoffset=1.0
!P.font = 0 

TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
!P.charsize = 1.5
!P.thick = 3
!P.charthick = 1

restore, 'mmintemplate.sav'
mmin_data = read_ascii('set22/mmin.out', template=mmintemplate) 
z = mmin_data.field1
mminc = mmin_data.field2*1.0e10
mminh = mmin_data.field3*1.0e10
mminav = mmin_data.field4*1.0e10

plot, z, mminav, /xlog, /ylog, xrange=[2,100], xtitle='redshift', $
      ytitle='M!Dmin!N (Msun)', xstyle=1, xthick=3, ythick=3 
oplot, z, mminc, linestyle=2, color=3, thick=3
oplot, z, mminh, linestyle=2, color=2, thick=3 

legend, ['H I', 'H II'], linestyle=[2,2], color=[3,2], $
        /bottom, box=0

device, /close_file
set_plot, 'X'

