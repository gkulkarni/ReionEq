# File: Makefile
# Cre: 2010-04-07
# Mod: $Date: 2012/09/17 07:25:40 $; $Revision: 1.25 $
#
# Compile reionization programs reion, qmass, pdelta 
# and massfn. Set FFLAGS to -g for debugging.

# FC =  /opt/intel/fce/10.1.018/bin/ifort
FC = gfortran 
FFLAGS = -g 
OBJECTS = constants.o storage.o interfaces.o srcov.o interpolate.o dqage.o \
	luminosity.o dqk61.o dqk51.o dqk41.o dqk31.o dqk21.o dqk15.o \
	counter.o dqagie.o nuint.o d1mach.o hub.o gammaph.o pigmd.o rtbis.o \
	zbrac.o getjmc.o getjmh.o getff.o lumfn.o gallum.o nfchk.o ztdyn.o \
	fpop3.o nuint_pop3.o nuint_pop2.o accrate.o getsfr.o outflow.o \
	interpolate2.o ejrate.o bi_interpolate2.o getmet.o haloyield.o\
	sfr_rollinde.o rtnewt.o calculate_nuintegral.o ejfrac.o outfrac.o\
	haloyield_species.o ngammafrac.o ejfrac_nonira.o outfrac_nonira.o\
	haloyield_nonira.o haloyield_species_nonira.o hallum.o hallum2.o scale_sed.o

REOBJ = constants.o storage.o dqage.o d1mach.o dqk61.o dqk51.o \
	dqk41.o dqk31.o dqk21.o dqk15.o dqpsrt.o

FORDELTA = storage.o constants.o interfaces.o counter.o zeroin.o \
	d1mach.o tabz.o

FORWL = constants.o storage.o interfaces.o srcwl.o interpolate.o dqage.o \
	luminosity.o dqk61.o dqk51.o dqk41.o dqk31.o dqk21.o dqk15.o \
	counter.o dqagie.o nuint.o d1mach.o hub.o gammaph.o pigmd.o rtbis.o \
	zbrac.o getmlow.o getjmc.o getjmh.o getff.o lumfn.o gallum.o 

FORHE = $(OBJECTS) 

FORSAMPLE = constants.o 

FORDELTAT = constants.o interfaces.o hub.o 

#readin.inc volav.inc 

%.o: %.f90
	$(COMPILE.f) $(OUTPUT_OPTION) $<

reion : $(OBJECTS) main.o readin.inc 
	$(FC) $(FFLAGS) -o reion $(OBJECTS) main.o

testejrate: $(OBJECTS) testejrate.o 
	$(FC) $(FFLAGS) -o testejrate $(OBJECTS) testejrate.o

re : $(REOBJ) re.o 
	$(FC) $(FFLAGS) -o re $(REOBJ) re.o 

pdelta : $(FORDELTA) pdelta.o 
	$(FC) $(CFLAGS) -o pdelta $(FORDELTA) pdelta.o

wl : $(FORWL) reion-wl.o 
	$(FC) $(CFLAGS) -o wl $(FORWL) reion-wl.o 

he : $(FORHE) halo_evolve.o 
	$(FC) $(CFLAGS) -o he $(FORHE) halo_evolve.o 

yld : $(OBJECTS) main_yld.o readin.inc 
	$(FC) $(FFLAGS) -o yld $(OBJECTS) main_yld.o

sample : $(FORSAMPLE) sample_ns.o 
	$(FC) $(FFLAGS) -o sample_ns $(FORSAMPLE) sample_ns.o 

deltat : $(FORDELTAT) deltat.o 
	$(FC) $(FFLAGS) -o deltat $(FORDELTAT) deltat.o 

.PHONY : clean
clean :
	-rm -rf *.mod *.o 

