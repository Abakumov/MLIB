# Makefile fuer Heidimod
# Erik H. Saenger
# 15.07.2008


# MPIRUN
MPIRUN=mpirun

# Fortran90 Compiler
F90=mpif90

# Fortran90 Optimierunsparameter
F90OPT = -O3 -msse3 -mcmodel=large 

# Fortran90 Endungen
F90EXT=f90
# Check all

go : moduli
	$(MPIRUN) -np 81 --hostfile ./machinefile  Evamod2014
	# STx zsnapheader &
	# STx xsnapheader &

# Clean

superclean :
	make clean
	rm *super*

clean :
	rm -f *.$(F90EXT)
	rm -f junk*
	rm -f zsnap*
	rm -f xsnap*
	rm -f ysnap*
	rm -f zseis*
	rm -f xseis*
	rm -f yseis*
	rm -f arb*
	rm -f geoloc
	rm -f PI*
	rm -f *~
	rm -f *.o
	rm -f coeff*
	rm -f wavelet*
	rm -f parfile
	rm -f *moduli*
	rm -f *header*
	rm -f slocation
	rm -f omegavisc yvisc*
	rm -f *super*
	rm -f me*
	rm -f trmfield
	rm -f sefield
	rm -f pefield
	rm -f gefield
	rm -f gefieldmean
	rm -f pefieldmean
	rm -f trmfieldmean
	rm -f sefieldmean



# This is only to create a model #+wavelet+coeff


moduli : modelmaker.F90
	cp modelmaker.F90 modelmaker.$(F90EXT)
	$(F90) -c $(F90OPT) modelmaker.$(F90EXT)
	$(F90) modelmaker.o -o modelmaker
	./modelmaker
	rm ./modelmaker

