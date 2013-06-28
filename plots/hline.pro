;+
; NAME:
;      HLINE
;     
; PURPOSE:
;      Draw a horizontal line on a pre-existing plot window.
;
; CALLING SEQUENCE:
;      HLINE, VAL
;
; INPUTS:
;
;      VAL: The y-value where the vertical line should be drawn
;
; KEYWORD PARAMETERS:
;
;      All keyword parameters are passed to OPLOT.
;
; SIDE EFFECTS:
;
;      Causes a horizontal line to appear on your screen.
;
; RESTRICTIONS:
;
;      This program won't do anything else. Sorry, them's the 
;      restrictions.
;
; EXAMPLE:
;
;      Draw a horizontal line at x = 4
;      IDL> plot, findgen(10)
;      IDL> hline, 4
;
; MODIFICATION HISTORY:
; Written sometime in 2003 by JohnJohn
; 4 June 2009: Modified to handle case of /xlog (hat tip to Michael Williams)
;
;-

pro hline, val,_extra=extra
if !x.type eq 1 then xrange = 10^!x.crange else xrange = !x.crange
nv = n_elements(val)
for i = 0, nv-1 do oplot,xrange,fltarr(2)+val[i],_extra=extra 
end


