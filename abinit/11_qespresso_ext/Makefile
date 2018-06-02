include ../make.inc

LIBNAME = lib11_qespresso_ext.a

SRC = \
flib_pwscf.F90  functionals_pwscf.F90  interfaces_11_qespresso_ext.F90 \
m_pseudo_pwscf.F90  read_upf_pwscf.F90

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


flib_pwscf.o : flib_pwscf.F90 
functionals_pwscf.o : functionals_pwscf.F90 flib_pwscf.o 
interfaces_11_qespresso_ext.o : interfaces_11_qespresso_ext.F90 
m_pseudo_pwscf.o : m_pseudo_pwscf.F90 
read_upf_pwscf.o : read_upf_pwscf.F90 interfaces_11_qespresso_ext.o m_pseudo_pwscf.o flib_pwscf.o


