include ../make.inc

LIBNAME = lib53_spacepar.a

SRC = \
dotprod_vn.F90  interfaces_53_spacepar.F90  meanvalue_g.F90  multipoles_fftr.F90


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


dotprod_vn.o : dotprod_vn.F90 
interfaces_53_spacepar.o : interfaces_53_spacepar.F90 
meanvalue_g.o : meanvalue_g.F90 
multipoles_fftr.o : multipoles_fftr.F90 
