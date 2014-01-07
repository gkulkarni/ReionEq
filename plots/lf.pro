FUNCTION schechter, z 

  ; Returns best-fit Schechter-fn. parameters.  
  ; See Kuhlen and Faucher-Giguere 2012 (table 1). 

  a_mstar = -20.42
  b_mstar = 0.27
  mstar = a_mstar + b_mstar*(z-6.0)

  a_phistar = -3.01
  b_phistar = -0.07
  phistar = a_phistar + b_phistar*(z-6.0)

  a_alpha = -1.84
  b_alpha = -0.06
  alpha = a_alpha + b_alpha*(z-6.0)

  param_array = fltarr(3) 
  param_array[0] = mstar
  param_array[1] = phistar
  param_array[2] = alpha

  return, param_array 
  
END 

PRO lf, set_name, row, opt 

  ; Plots luminosity function from our Pop3 model in variety of ways.
  window, xsize=750, ysize=750
  Device, decomposed=0
  TvLCT, 255, 0, 0, 2 
  !P.charsize = 2.0

  z_init = 49.5
  dz = -0.5 
  z = row*dz + z_init 

  restore, 'halo_template.sav'
  n_data = read_ascii(strtrim(set_name)+'/nofmc.out', template=stars_template)
  m_data = read_ascii(strtrim(set_name)+'/halos.out', template=stars_template)
  l_data = read_ascii(strtrim(set_name)+'/l1500.out', template=stars_template)
  mstrdot_data = read_ascii(strtrim(set_name)+'/halos_aux.out', template=stars_template)

  n = n_data.field001[*,row]
  m = m_data.field001[*,row]*1.0d10
  l = l_data.field001[*,row]

  total_indices = size(n, /n_elements)
  index_lowerlimit = 0 
  for i = 1, total_indices do begin 
     if (l[i-1] lt -1.0) then begin 
        index_lowerlimit = i-1
        break 
     endif
  endfor

  index_upperlimit = total_indices-1 
  for i = index_lowerlimit+1, total_indices do begin 
     if (l[i-1] gt -1.0) then begin 
        index_upperlimit = i-2
        break 
     endif
  endfor

  n = n_data.field001[index_lowerlimit:index_upperlimit,row]
  m = m_data.field001[index_lowerlimit:index_upperlimit,row]*1.0d10
  l = l_data.field001[index_lowerlimit:index_upperlimit,row]
  mstrdot = mstrdot_data.field001[index_lowerlimit:index_upperlimit,row]*1.0e10 
  l = smooth(l, 80, /edge_truncate)

  n_indices = size(n, /n_elements)
  schechter_params = schechter(z) 
  mstar = schechter_params[0]
  phistar = 10.0^schechter_params[1]
  alpha = schechter_params[2]
  lf_schechter = fltarr(n_indices)
  lf = fltarr(n_indices) 
  for i = 0, n_indices-1 do begin 
     lf_schechter[i] = 0.4*alog(10.0)*phistar*((10.0^(0.4*(mstar-l[i])))^(alpha+1))*$
                       exp(-10.0^(0.4*(mstar-l[i]))) 

     if (i eq n_indices-1) then begin 
        dl = 1.0                ; just a test value 
     endif else begin 
                                ; Divide by dl to change from mass function to luminosity fn.
        dl = abs(l[i+1]-l[i]) 
        lf[i] = n[i]/dl 
     endelse
  endfor

                                ; Below I create a toy luminosity assignment to study the LF. 
                                ; Try opt=4 to understand the difference from the model prediction.
  l_fit = fltarr(n_indices)
  lf_fit = fltarr(n_indices)
  lf_schechter_fit = fltarr(n_indices)
  for i = 0, n_indices-1 do begin 
                                ; This is a linear fit. 
     l_fit[i] = ((l[10]-l[12])/(10.0-12.0))*(float(i)-10.0) + l[12]
     lf_schechter_fit[i] = 0.4*alog(10.0)*phistar*(10.0^(0.4*(mstar-l_fit[i])*(alpha+1)))*$
                           exp(-10.0^(0.4*(mstar-l_fit[i]))) 
     if (i eq n_indices-1) then begin 
        dl = 1.0                ; just a test value 
     endif else begin 
        dl = abs(l_fit[i+1]-l_fit[i]) 
        lf_fit[i] = n[i]/dl 
     endelse
  endfor

  ; readcol, '../result/phi_tau0.09.dat', z_2010, mag_2010, lf_2010, /silent 
  readcol, '../result/lumfn_z7.0.dat', z_2010, mag_2010, lf_2010, /silent 
  readcol, 'nofm_2010.chk', m_2010, nm_2010, nmalt_2010, /silent 

  readcol, 'phi_z7.0.dat', mab_z7, phi_z7, philo_z7, phihi_z7 
  phi_z7 = 10.0^phi_z7 

  n_indices = size(nm_2010,/n_elements)
  for i = 0, n_indices-1 do begin 
     if (i eq n_indices-1) then begin 
        dm = 1.0                ; just a test value 
     endif else begin 
        dm = m_2010[i+1]-m_2010[i]
     endelse

     nm_2010[i] = nm_2010[i]*dm       ; Mpc^-3
     nmalt_2010[i] = nmalt_2010[i]*dm ; Mpc^-3
  endfor
  m_2010 = m_2010 * 1.0e10      ; Msun 

  case opt of 
     1: plot, m, n, /xlog, /ylog, xtitle="halo mass (M!D!9n!X!N)", ytitle="N (Mpc!E-3!N)"
     2: plot, m, l, /xlog, xtitle="halo mass (M!D!9n!X!N)", ytitle="M!D1500!N"
     3: plot, l, n, /ylog, xtitle="M!D1500!N", ytitle="N (Mpc!E-3!N)"
     4: begin
        plot, l, xtitle="array index", ytitle="M!D1500!N"
        oplot, l_fit, linestyle=2 
     end
     5: begin 
        plot, l, n, /ylog, xtitle="M!D1500!N", ytitle="N (Mpc!E-3!N)"
        oplot, l, lf_schechter, linestyle=2 
     end
     6: plot, mstrdot 
     7: plot, l, mstrdot, /ylog 
     8: begin 
        plot, l, lf, /ylog, yrange=[1.0e-6,10.0]
        oplot, l, lf_schechter, linestyle=2 
        oplot, l, lf_fit, color=2
     end
     9: plot, n, /ylog
     10: plot, l_fit, lf_fit, /ylog, yrange=[1.0e-6,10.0]
     11: begin
        plot, l_fit, lf_fit, /ylog, xrange=[-20.0,-12.0], yrange=[1.0e-6,1.0]
        oplot, l_fit, lf_schechter_fit, linestyle=2
     end
     12: begin 
        plot, l_fit, lf_fit, /ylog, xrange=[-22.0,-12.0], yrange=[1.0e-6,10.0], $
              xtitle='!6M!D1500!N', ytitle='LF' 
        oplot, l_fit, lf_schechter_fit, linestyle=2
        oplot, mag_2010, lf_2010, linestyle=1 
        oplot, l, lf, color=2
        oplot, mab_z7, phi_z7, psym=6, symsize=2 
        legend, ['model', 'data', '2010 result', 'lin. fit'], $
                linestyle=[0,2,1,0], color=[2,-1,-1,-1], /bottom, /right
     end
     13: begin 
        plot, m_2010, nm_2010, /xlog, /ylog  
        oplot, m, n, color=2, thick=4
     end
     14: begin 
        plot, l_fit, lf_fit, /ylog, xrange=[-25.0,-12.0], yrange=[1.0e-6,1.0], $
              xtitle='!6M!D1500!N', ytitle='LF', /nodata
        oplot, l_fit, lf_schechter_fit, linestyle=2
        oplot, l, lf, color=2
        legend, ['model', 'data'], linestyle=[0,2], color=[2,-1], /bottom, /right
     end
  endcase

END 

