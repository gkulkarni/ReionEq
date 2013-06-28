; Plot our result. 
restore, 'sfrfiletemplate.sav'
sfrdata = read_ascii('set76/sfr.out', template=sfrfiletemplate)
sfr_tot = sfrdata.total_sfr ; Msun yr^-1 Mpc^-3 
z = sfrdata.redshift
plot, z, sfr_tot, /xlog, /ylog, xrange=[1,100]

sfrdata = read_ascii('set82/sfr.out', template=sfrfiletemplate)
sfr_tot = sfrdata.total_sfr ; Msun yr^-1 Mpc^-3 
z = sfrdata.redshift
oplot, z, sfr_tot, color=4
