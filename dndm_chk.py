#!/usr/bin/python 

# File: dndm_chk.py
#  Cre: 2012
#  Mod: $Date$ ($Revision$)

# This small program integrates the SEDs over frequency ranges to
# estimate the total fluxes from stars, just as a check.

import sys 
sedfile = open(sys.argv[1],"r")
record = sedfile.readline()
(nu, dndmdnu) = record.split()
nu1 = float(nu) 
dndmdnu1 = float(dndmdnu) 
dndm = 0.0 
for record in sedfile: 
    (nu, dndmdnu) = record.split()
    nu2 = float(nu) 
    dndmdnu2 = float(dndmdnu) 
    dndm = dndm + 0.5*(nu1-nu2)*(dndmdnu1+dndmdnu2) 
    dndmdnu1 = dndmdnu2 
    nu1 = nu2 
print dndm 
sedfile.close()
    
