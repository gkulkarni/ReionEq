
; File: sfr_ratio.f90
;  Cre: 12-2012
;  Mod: $Date: 2013/01/18 17:10:23 $ ($Revision: 1.1 $) 

; Plots the Pop II and Pop III SFR from reion output.  Also plots
; their ratio and observational data points. 

;; set_plot, 'ps'
;; device, filename='sfr_ratio.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

window, xsize=1000, ysize=1000
Device, decomposed=0
!P.multi = [0,1,2,0,1]
TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 0, 255, 4
!P.charsize = 1.2
!P.charthick = 1
!P.thick = 1 

readcol, '../data/sfrsfr.data', rs, rserr, sf, sferr, /silent  
sfe = 10.0^sf 
plotsym, 0, 0.5, /FILL
plot, rs, sfe, psym=8, /ylog, /xlog, xrange=[5.0e-1,100], $
      ytitle='SFR (M!D!9n!N!X yr!E-1!N Mpc!E-3 !N)', xstyle=1, yrange=[1.0e-10,1], $
      ytickformat='Exp1', position=[0.1,0.1,0.9,0.7], xtitle='!6z', /nodata 
oplot, rs, sfe, psym=8, color=2
yup = 10.0^(sf+sferr*0.5)
ylo = 10.0^(sf-sferr*0.5)
dy = yup-ylo
plot_err, rs, sfe, dy, dx1=rserr, color=2

restore, 'sfrfiletemplate.sav'
sfrdata = read_ascii('set55/sfr.out', template=sfrfiletemplate)
sfr_tot = sfrdata.total_sfr
z = sfrdata.redshift
;oplot, z, sfr_tot 
sfr_pop2 = sfrdata.pop2_sfr
sfr_pop3 = sfrdata.pop3_sfr
n = size(sfr_pop3, /n_elements)
sfr_pop3[84]=1.0e-15
oplot, z, sfr_pop3, linestyle=5
pop3_frac1 = sfr_pop3/sfr_tot

sfrdata = read_ascii('set7/sfr.out', template=sfrfiletemplate)
sfr_tot2 = sfrdata.total_sfr
oplot, z, sfr_tot2, color=2
oplot, z, sfr_tot
z = sfrdata.redshift
sfr_pop2 = sfrdata.pop2_sfr
sfr_pop3 = sfrdata.pop3_sfr
sfr_pop3[84]=1.0e-15
oplot, z, sfr_pop3, linestyle=5, color=2
pop3_frac2 = sfr_pop3/sfr_tot2

vline, 8.5, linestyle=2
xyouts, 6.2, 1.0e-3, 'z!Dreion!N', orientation=90.0, charsize=1.5, alignment=0.5

; Add legend.
legend, ['Hopkins+Beacom 06','1-100 M!D!9n!X!N (total)', '1-100 M!D!9n!X!N (pop. III)', $
         '100-260 M!D!9n!X!N (total)', '100-260 M!D!9n!X!N (pop. III)'], linestyle=[0,0,5,0,5], $
        psym=[8,0,0,0,0], color=[2,-1,-1,2,2], /bottom, charsize=1

;-----------------------------------------------------------------

tick = replicate(' ',3)
plot, z, pop3_frac1, /ylog, /xlog, xrange=[5.0e-1,100], xstyle=1, yrange=[1.0e-4,1.0], $
      ytickformat='Exponent', position=[0.1,0.7,0.9,0.97], xtickname=tick, $
      ytitle='ratio (popIII/total)'
oplot, z, pop3_frac2, color=2

vline, 8.5, linestyle=2
xyouts, 6.2, 1.0e-2, 'z!Dreion!N', orientation=90.0, charsize=1.5, alignment=0.5

; Add legend.
legend, ['1-100 M!D!9n!X', '100-260 M!D!9n!X'], linestyle=[0,0], $
        color=[-1,2], /bottom, charsize=1

;-----------------------------------------------------------------

;; device, /close_file
;; set_plot, 'X'

