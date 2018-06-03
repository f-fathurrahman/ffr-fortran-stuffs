include ../make.inc

LIBNAME = lib62_cg_noabirule.a

SRC = \
bracketing.F90  cgpr.F90        dummy_functions.inc  frskerker2.F90  linmin.F90 \
brent.F90       dotproduct.F90  frskerker1.F90       interfaces_62_cg_noabirule.F90


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



bracketing.o : bracketing.F90 dummy_functions.inc 
brent.o : brent.F90 dummy_functions.inc interfaces_62_cg_noabirule.o 
cgpr.o : cgpr.F90 dummy_functions.inc interfaces_62_cg_noabirule.o 
dotproduct.o : dotproduct.F90 
frskerker1.o : frskerker1.F90 interfaces_62_cg_noabirule.o 
frskerker2.o : frskerker2.F90 interfaces_62_cg_noabirule.o 
interfaces_62_cg_noabirule.o : interfaces_62_cg_noabirule.F90 
linmin.o : linmin.F90 dummy_functions.inc interfaces_62_cg_noabirule.o 
