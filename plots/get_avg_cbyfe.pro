
PRO getav_cbyfe, set 

z = 2.0 
dz = 0.5 

repeat begin
   row = uint((50.0-z)/dz)-1 
   abr3, set, row, z 
   z = z + 0.5 
endrep until z gt 10.5

END

