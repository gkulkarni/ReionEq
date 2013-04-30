#!/usr/bin/python

import sys
ratesfile = open(sys.argv[1],"r")
sfrfile = open(sys.argv[2],"r")
sfrdata = sfrfile.readlines()
i = 0 
for record in ratesfile:
    entry = record.split()
    ab = entry[1]
    o = entry[2]
    e = entry[3]
    sfrentry = sfrdata[i].split()
    sfr = sfrentry[3]
    i = i+1 
    print float(entry[0]), float(ab)+float(e), float(o)+(float(sfr)*1.0e-10)
    
    
    
