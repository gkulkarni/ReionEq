function metlabel, axis, index, number 

  if number eq 0 then begin 
     return, ' ' 
  endif else begin 
     return, string(number, format='(I2)')
  endelse

end




