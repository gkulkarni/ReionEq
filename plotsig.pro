sig = fltarr(2,100)
openr, lun, '../foo', /get_lun
readf, lun, sig
d = sig[0,*]
p = sig[1,*]
plothist, d, /autobin
