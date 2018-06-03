include ../make.inc

LIBNAME = lib78_eph.a

SRC = \
interfaces_78_eph.F90  m_gkk.F90  m_phgamma.F90  m_phpi.F90  m_sigmaph.F90


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



interfaces_78_eph.o : interfaces_78_eph.F90 
m_gkk.o : m_gkk.F90 
m_phgamma.o : m_phgamma.F90 
m_phpi.o : m_phpi.F90 
m_sigmaph.o : m_sigmaph.F90 
