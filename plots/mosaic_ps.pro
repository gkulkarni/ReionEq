PRO mosaic_ps, set_lowmass, set_highmass 

; set_lowmass: result set corresponding to 1-100 Msun IMF. 
; set_highmass: result set corresponding to 100-260 Msun IMF. 

set_plot, 'ps'
device, filename='mosaic.ps', xsize=24.0, ysize=8.0, /inches, color=1, yoffset=1.0

TvLCT, 255, 0, 0, 2 
TvLCT, 0, 127, 255, 3
TvLCT, 255, 255, 0, 4
TvLCT, 0, 255, 0, 5
TvLCT, 255, 255, 255, 6
!P.charsize = 3.0
!P.thick = 1.0
!P.multi = [0,6,1,0,0]

plotwidth = 0.17 
plotxspace = 0.06
x1 = 0.05 
x2 = x1 + plotwidth + plotxspace 
x3 = x2 + plotwidth + plotxspace 

; Plot 1: Single halo Z evolution. 

restore, 'mmintemplate.sav'
mmin_data = read_ascii(set_highmass + '/mmin.out', template=mmintemplate) 
z = mmin_data.field1
mminh = alog10(mmin_data.field3)+10.0
tick = replicate(' ',3)
plot, z, mminh, /xlog, xrange=[1,50], linestyle=5, ytitle = 'log!D10!N(M/M!D!9n!X!N)', $
      yrange=[6,12], position=[x1,0.7,x1+plotwidth,0.97], xtickname=tick, xstyle=1
restore, 'halo_template.sav'
halos_data = read_ascii(set_highmass + '/halos.out', template=stars_template)
halos = alog10(halos_data.field001[170,*])+10.0
oplot, z, halos, color=2
halos = alog10(halos_data.field001[140,*])+10.0
oplot, z, halos, color=3
halos = alog10(halos_data.field001[100,*])+10.0
oplot, z, halos, color=5

legend, ['!NM!Dmin!N (HII region)'], linestyle=[5], charsize=1.1, /bottom, number=3
legend, ['(a1)'], /right, charsize=1.1, box=0 

metal_data = read_ascii(set_highmass + '/halos_metals.out', template=stars_template)
gas_data =  read_ascii(set_highmass + '/halos_gas.out', template=stars_template)
metal = metal_data.field001[100,*]
gas = gas_data.field001[100,*]
Zmetal = fltarr(100)
for i = 0, 99 do begin 
   Zmetal[i] = (metal[0,i]/gas[0,i])/0.02
   print, z[i], Zmetal[i], metal[0,i], gas[0,i] 
endfor 
Zmetal = alog10(Zmetal) 
tick = replicate(' ',3)
plot, z, Zmetal, /xlog, xrange=[1,50], ytitle='log!D10!N(Z/Z!D!9n!X!N)', $
      xstyle=1, yrange=[-6.0,0.0], /nodata, position=[x1,0.1,x1+plotwidth,0.7], $
      xtitle='redshift', ytickformat='metlabel', xtickname=tick
axlabel, [1.0, 10.0, 50.0], /xaxis, charsize=1.5, format='(I)'
oplot, z, Zmetal, color=5 

legend, ['Pop. III IMF: 100-260 M!D!9n!X', 't!Denrich!N=10!E9!N Myr'], charsize=1.1, box=0
; legend, ['Pop. III IMF: 100-260 M!D!9n!X', 't!Denrich!N=0'], charsize=1.1, box=0
legend, ['(a2)'], /right, charsize=1.1, box=0 

metal_data = read_ascii(set_highmass + '/halos_metals.out', template=stars_template)
gas_data =  read_ascii(set_highmass + '/halos_gas.out', template=stars_template)
metal = metal_data.field001[140,*]
gas = gas_data.field001[140,*]
Zmetal = fltarr(100)
for i = 0, 99 do begin 
   Zmetal[i] = (metal[0,i]/gas[0,i])/0.02
endfor 
Zmetal = alog10(Zmetal) 
oplot, z, Zmetal, color=3

metal_data = read_ascii(set_highmass + '/halos_metals.out', template=stars_template)
gas_data =  read_ascii(set_highmass + '/halos_gas.out', template=stars_template)
metal = metal_data.field001[170,*]
gas = gas_data.field001[170,*]
Zmetal = fltarr(100)
for i = 0, 99 do begin 
   Zmetal[i] = (metal[0,i]/gas[0,i])/0.02
endfor 
Zmetal = alog10(Zmetal) 
oplot, z, Zmetal, color=2

hline, -4.0, linestyle=2
xyouts, 2.0, -3.8, 'Z!Dcrit!N', charsize=1.5, alignment=0.5

; Plot GPI 

restore, 'reionfiletemplate_splitgpi.sav'
reiondata = read_ascii(set_lowmass + '/reion.out', template=reionfiletemplate_splitgpi)
redshift = reiondata.field01
gpi = reiondata.field04
gpi2 = reiondata.field05
gpi3 = reiondata.field06
gpi = gpi * 1.0e12 
gpi2 = gpi2 * 1.0e12 
gpi3 = gpi3 * 1.0e12 
gpi3 = smooth(gpi3, 3) 

ratio1 = gpi3/gpi

tick = replicate(' ',3)
plot, redshift, gpi, /ylog, xrange=[1,50], yrange=[1.0e-10, 1.0e6], xstyle=1, xtitle='!6redshift', $
      ytitle='log!D10!N(!7C!6!DHI!N/10!E-12!Ns!E-1!N)', ytickformat='exp2',$
      /xlog, /nodata, position=[x2,0.1,x2+plotwidth,0.7], xtickname=tick, ystyle=1
axlabel, [1.0, 10.0, 50.0], /xaxis, charsize=1.5, format='(I)'
oplot, redshift, gpi3, linestyle=5 

reiondata = read_ascii(set_highmass + '/reion.out', template=reionfiletemplate_splitgpi)
redshift = reiondata.field01
gpi = reiondata.field04
gpi2 = reiondata.field05
gpi3 = reiondata.field06
gpi = gpi * 1.0e12 
gpi2 = gpi2 * 1.0e12 
gpi3 = gpi3 * 1.0e12 
gpi3 = smooth(gpi3, 3) 

ratio2 = gpi3/gpi 

oplot, redshift, gpi, color=2
oplot, redshift, gpi3, color=2, linestyle=5 

plotsym, 0, 0.5, /FILL

readcol, '../data/gammapi_mw.dat', x, y, dy1, dy2, /silent 
oplot, x, y, psym=8, color=3
oploterror, x, y, dy1, errcolor=3, psym=3, /hibar, /nohat
oploterror, x, y, dy2, errcolor=3, psym=3, /lobar, /nohat 
x0 = x[0]
y0 = y[0] 
x1 = x[0]
y1 = 7.5e-2

readcol, '../data/gammapi_bh.dat', x, y, dy1, dy2, /silent 
x[0] = 0.0
oplot, x, y, psym=8, color=2
oploterror, x, y, dy1, errcolor=2, psym=3, /hibar, /nohat
oploterror, x, y, dy2, errcolor=2, psym=3, /lobar, /nohat
x0 = x[4]
y0 = y[4] 
x1 = x[4]
y1 = 1.0e-1 
arrow, x0, y0, x1, y1, /data, hsize=7.0, color=2 

readcol, '../data/gammapi_cafg.dat', x, y, dy, /silent 
x[0] = 0.0 
oplot, x, y, psym=8, color=5
oploterror, x, y, dy, errcolor=5, psym=3, /nohat 

vline, 8.5, linestyle=2
xyouts, 6.5, 1.0e-8, 'z!Dreion!N', orientation=90.0, charsize=1.5, alignment=0.5

al_legend, ['Faucher-Giguere 08', 'Meiksin and White 04', 'Bolton and Haehnelt 07', $
            '1-100 M!D!9n!X!N (total)', '1-100 M!D!9n!X!N (pop. III)', '100-260 M!D!9n!X!N (total)', $
            '100-260 M!D!9n!X!N (pop. III)'], linestyle=[0,0,0,0,5,0,5], psym=[8,8,8,0,0,0,0], $
           color=[5,3,2,-1,-1,2,2], /right, charsize=1.1, background_color=6
legend, ['(b2)'], /bottom, charsize=1.1, box=0 

tick = replicate(' ',3)
plot, redshift, ratio1, position=[x2,0.7,x2+plotwidth,0.97], /xlog, /ylog, xrange=[1,50], $
      xstyle=1, xtickname=tick, ytitle='ratio', ytickformat='Exponent', yrange=[1.0e-4,1.0]
oplot, redshift, ratio2, color=2 
vline, 8.5, linestyle=2
xyouts, 6.5, 1.0e-2, 'z!Dreion!N', orientation=90.0, charsize=1.5, alignment=0.5
legend, ['(b1)'], charsize=1.1, box=0 
al_legend, ['1-100 M!D!9n!X', '100-260 M!D!9n!X'], linestyle=[0,0], $
        color=[-1,2], /bottom, charsize=1.1, background_color=6


; Plot 3: SFR 

readcol, '../data/sfrsfr.data', rs, rserr, sf, sferr, /silent  
sfe = 10.0^sf 
plotsym, 0, 0.5, /FILL
tick = replicate(' ',3)
plot, rs, sfe, psym=8, /ylog, /xlog, xrange=[1,50], $
      ytitle='SFRD (M!D!9n!N!X yr!E-1!N Mpc!E-3 !N)', xstyle=1, yrange=[1.0e-10,1], $
      ytickformat='Exp1', position=[x3,0.1,x3+plotwidth,0.7], xtitle='!6redshift', $
      /nodata, xtickname=tick
axlabel, [1.0, 10.0, 50.0], /xaxis, charsize=1.5, format='(I)'
oplot, rs, sfe, psym=8, color=2
yup = 10.0^(sf+sferr*0.5)
ylo = 10.0^(sf-sferr*0.5)
dy = yup-ylo
plot_err, rs, sfe, dy, dx1=rserr, color=2

restore, 'sfrfiletemplate.sav'
sfrdata = read_ascii(set_lowmass + '/sfr.out', template=sfrfiletemplate)
sfr_tot = sfrdata.total_sfr
z = sfrdata.redshift
sfr_pop2 = sfrdata.pop2_sfr
sfr_pop3 = sfrdata.pop3_sfr
n = size(sfr_pop3, /n_elements)
sfr_pop3[84]=1.0e-15
oplot, z, sfr_pop3, linestyle=5
pop3_frac1 = sfr_pop3/sfr_tot

sfrdata = read_ascii(set_highmass + '/sfr.out', template=sfrfiletemplate)
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

al_legend, ['Hopkins+Beacom 06','1-100 M!D!9n!X!N (total)', '1-100 M!D!9n!X!N (pop. III)', $
         '100-260 M!D!9n!X!N (total)', '100-260 M!D!9n!X!N (pop. III)'], linestyle=[0,0,5,0,5], $
        psym=[8,0,0,0,0], color=[2,-1,-1,2,2], /bottom, charsize=1.1, background_color=6
legend, ['(c2)'], /right, charsize=1.1, box=0 

tick = replicate(' ',3)
plot, z, pop3_frac1, /ylog, /xlog, xrange=[1,50], xstyle=1, yrange=[1.0e-4,1.0], $
      ytickformat='Exponent', position=[x3,0.7,x3+plotwidth,0.97], xtickname=tick, $
      ytitle='ratio (popIII/total)'
oplot, z, pop3_frac2, color=2

vline, 8.5, linestyle=2
xyouts, 6.2, 1.0e-2, 'z!Dreion!N', orientation=90.0, charsize=1.5, alignment=0.5

al_legend, ['1-100 M!D!9n!X', '100-260 M!D!9n!X'], linestyle=[0,0], $
        color=[-1,2], /bottom, charsize=1.1, background_color=6
legend, ['(c1)'], charsize=1.1, box=0 

device, /close_file
set_plot, 'X'

end

