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

  for i = index_lowerlimit+1, total_indices do begin 
     if (l[i-1] gt -1.0) then begin 
        index_upperlimit = i-2
        break 
     endif
  endfor

  n = n_data.field001[index_lowerlimit:index_upperlimit,row2]
  m = m_data.field001[index_lowerlimit:index_upperlimit,row2]*1.0d10
  l = l_data.field001[index_lowerlimit:index_upperlimit,row2]
  l = smooth(l, 80, /edge_truncate)

  case opt of 
     1: oplot, m, n, color=2
     2: oplot, m, l, color=2
     3: oplot, l, n, color=2
     4: oplot, l, color=2
  endcase

