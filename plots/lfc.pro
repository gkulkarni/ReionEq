PRO lfc, set_name, row, row2, opt, incr 

  window, xsize=1000, ysize=1000
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
  mstrdot_data = read_ascii(strtrim(set_name)+'/halos_mstardot.out', template=stars_template)
  sfrcontrib_data = read_ascii(strtrim(set_name)+'/halosh_sfrcontrib.out', template=stars_template)

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

  print, 'index_lowerlimit = ', index_lowerlimit 

  index_upperlimit = total_indices-1 
  for i = index_lowerlimit+1, total_indices do begin 
     if (l[i-1] gt -1.0) then begin 
        index_upperlimit = i-2
        break 
     endif
  endfor

  n = n_data.field001[index_lowerlimit:index_upperlimit,row]
  m = m_data.field001[index_lowerlimit:index_upperlimit,row]*1.0d10
  l = l_data.field001[index_lowerlimit:index_upperlimit,row]+incr
  mstrdot = mstrdot_data.field001[index_lowerlimit:index_upperlimit,row]*1.0e10 
  sfrcontrib = sfrcontrib_data.field001[index_lowerlimit:index_upperlimit,row]*1.0e10 

  case opt of 
     1: plot, m, n, /xlog, /ylog, xtitle="halo mass (M!D!9n!X!N)", ytitle="N (Mpc!E-3!N)"
     2: plot, m, l, /xlog, xtitle="halo mass (M!D!9n!X!N)", ytitle="M!D1500!N"
     3: plot, l, n, /ylog, xtitle="M!D1500!N", ytitle="N (Mpc!E-3!N)"
     4: plot, l, xtitle="array index", ytitle="M!D1500!N"
     5: begin 
        l = smooth(l, 80, /edge_truncate)
        plot, l, n, /ylog, xtitle="M!D1500!N", ytitle="N (Mpc!E-3!N)"
     end
     6: plot, mstrdot 
     7: plot, l, mstrdot, /ylog 
     8: plot, l, sfrcontrib, /ylog 
     9: plot, m, /ylog
  endcase

  n = n_data.field001[*,row2]
  m = m_data.field001[*,row2]*1.0d10
  l = l_data.field001[*,row2]

  total_indices = size(n, /n_elements)

  for i = 1, total_indices do begin 
     if (l[i-1] lt -1.0) then begin 
        index_lowerlimit = i-1
        break 
     endif
  endfor

  print, 'index_lowerlimit = ', index_lowerlimit 

  for i = index_lowerlimit+1, total_indices do begin 
     if (l[i-1] gt -1.0) then begin 
        index_upperlimit = i-2
        break 
     endif
  endfor

  n = n_data.field001[index_lowerlimit:index_upperlimit,row2]
  m = m_data.field001[index_lowerlimit:index_upperlimit,row2]*1.0d10
  l = l_data.field001[index_lowerlimit:index_upperlimit,row2]+incr
  mstrdot = mstrdot_data.field001[index_lowerlimit:index_upperlimit,row2]*1.0e10 
  sfrcontrib = sfrcontrib_data.field001[index_lowerlimit:index_upperlimit,row]*1.0e10 

  case opt of 
     1: oplot, m, n, color=2
     2: oplot, m, l, color=2
     3: oplot, l, n, color=2
     4: oplot, l, color=2
     5: begin 
        l = smooth(l, 80, /edge_truncate)
        oplot, l, n, color=2
     end
     6: oplot, mstrdot, color=2 
     7: oplot, l, mstrdot, color=2
     8: oplot, l, sfrcontrib, color=2
     9: oplot, m, color=2
  endcase

END 

