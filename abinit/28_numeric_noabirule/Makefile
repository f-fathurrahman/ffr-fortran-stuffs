include ../make.inc

LIBNAME = lib28_numeric_noabirule.a

SRC = \
interfaces_28_numeric_noabirule.F90 \
m_numeric_tools.F90 \
m_slk.F90 \
cspint.F90  \
j0y0.F90      m_abilasi.F90     m_elpa.F90 \
dzgedi.F90  j1y1.F90      m_abi_linalg.F90  m_gaussian_quadrature.F90  m_sort.F90 \
dzgefa.F90  jacobi.F90    m_array.F90       m_lebedev.F90              m_splines.F90 \
gamma_function.F90  k0.F90  m_blas.F90  m_levenberg_marquardt.F90 \
k1.F90  m_bspline.F90  uniformrandom.F90 \
interp.F90  lib_four.F90  m_cgtools.F90     m_per_cond.F90 \
intrpl.F90  ludcmp.F90    m_cplxtools.F90   m_random_zbq.F90

# test_spline_integrate.F90

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



