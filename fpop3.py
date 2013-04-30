#!/usr/bin/python

# Cre: 12-2011
# Mod: $Date$ ($Revision$) 

# Read PINOCCHIO output and calculate the fraction of halos and given
# redshift that host Pop III stars.  See Schneider 2006.  Takes two
# command line arguments.

import sys 
import math 
if (len(sys.argv) != 3):
    print("Usage: ./fpop3.py <run_index> <redshift>")
    print(" E.g.: ./fpop3.py lcdm05 23.0") 

# Read PINOCCHIO merger tree data from frag.*.*.halo.out.
if (int(sys.argv[2]) > 9) : 
    filename = "./" + sys.argv[1] + "/" + "frag." + sys.argv[2] + ".0." + sys.argv[1] + ".halo.out"
else : 
    filename = "./" + sys.argv[1] + "/" + "frag." + sys.argv[2] + ".00." + sys.argv[1] + ".halo.out"
datafile = open(filename, "r")
data = datafile.readlines()

# Store data in a dictionary.
halos = {} 
for i in range(len(data)): 
    one_halo = data[i].rsplit()
    halos[one_halo[0]] = one_halo[1:]

particle_mass = 5.103861e6 # M_solar 

# Function `is_pop3' returns TRUE if halo with id `halo_id' is a Pop
# III halo and FALSE otherwise.
def check_pop3(halo_id): 
    halo_data = halos[halo_id]
    if (halo_data[1]==halo_id or halo_data[1]==candidate): 
        is_pop3 = True 
        return(is_pop3) 
    progenitor_data = halos[halo_data[1]]
    progenitor_mass = float(progenitor_data[3])*particle_mass 
    zprog = float(progenitor_data[len(progenitor_data)-1])
    mmin = 1.0e8*(1+zprog/10.0)**(-1.5)
    if (progenitor_mass > mmin): 
        is_pop3 = False 
        return(is_pop3) 
    else:
        return(check_pop3(halo_data[1]))

# Count haloes that have formed in the current redshift window.
new_halos = 0 
pop3_halos = 0 
for k, v in halos.iteritems():
    if (float(v[6]) == -1.0 and
        ((float(v[7]) == -1.0 and v[8] < (float(sys.argv[2])+1.0)) or
         (float(v[7]) != -1.0 and v[7] < (float(sys.argv[2])+1.0)))):
        new_halos = new_halos + 1 
        candidate = k 
        if (check_pop3(candidate)): 
            pop3_halos = pop3_halos + 1 

nbin_perdecade = 10 
m_high = 1.0e12
m_low = 1.0e7
nbin = int(math.log10(m_high/m_low)*nbin_perdecade)
ldm = 10.0**(1.0/float(nbin_perdecade))

mass_bins = {}
for i in range(nbin-1): 
    m2 = m_low * (ldm**i) 
    m1 = m_low * (ldm**(i-1))
    mass_bins[i] = [math.sqrt(m1*m2),0]
    pop3_halos = 0 
    for k, v in halos.iteritems():
        halo_mass = float(v[3])*particle_mass 
        halo_zform = float(v[len(v)-1])
        if ((halo_mass>=m1 and halo_mass<m2) and (halo_zform < (float(sys.argv[2])+1.0))): 
            mass_bins[i][1] = mass_bins[i][1] + 1 
            candidate = k 
            if (check_pop3(candidate)): 
                pop3_halos = pop3_halos + 1 
    if(mass_bins[i][1]!=0): 
        fp3 = float(pop3_halos)/float(mass_bins[i][1]) 
        print mass_bins[i][0], fp3
    else:
        print mass_bins[i][0], 0.0

