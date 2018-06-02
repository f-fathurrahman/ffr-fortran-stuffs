include ../make.inc

LIBNAME = lib62_poisson.a

SRC = \
interfaces_62_poisson.F90  psolver_hartree.F90  psolver_kernel.F90  psolver_rhohxc.F90


OBJ = $(SRC:.F90=.o)

#
# Suffix rule for Fortran 90
#
%.mod :
	@if [! -f $@ ]; then \
		rm $(*F).o; \
		fi
	$(MAKE) $<

%.o : %.F90
	$(F90) $(F90_OPTS) -c -o $(*F).o $<

#
# Fortran 77 sources
# supress warning
.SUFFIXES: .o .f
.f.o:
	$(F77) -c $(F77_OPTS) $<

#
# C source
#
.SUFFIXES: .o .c
.c.o:
	$(CC) -c $(CC_OPTS) $<

# Targets
lib: $(OBJ)
	ar rcs $(LIBNAME) *.o; mv $(LIBNAME) ../libs; mv *.mod ../mods

# does not delete *.x files
clean:
	rm -rf *.o *.mod $(LIBNAME)



interfaces_62_poisson.o : interfaces_62_poisson.F90 
psolver_hartree.o : psolver_hartree.F90 interfaces_62_poisson.o 
psolver_kernel.o : psolver_kernel.F90 
psolver_rhohxc.o : psolver_rhohxc.F90 interfaces_62_poisson.o 
