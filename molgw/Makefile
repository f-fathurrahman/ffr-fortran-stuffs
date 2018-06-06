include make.inc

LIBS_EXT =


SRC = \
lebedev_quadrature.f \
m_definitions.f90 \
m_gos.f90 \
m_warning.f90 \
m_tools.f90 \
m_mpi_auxil.f90  m_mpi_grid.f90  m_mpi_local.f90  m_mpi_ortho.f90  m_mpi_trans.f90  m_mpi_world.f90 m_mpi.f90 \
m_cart_to_pure.f90 \
m_gaussian.f90 \
m_elements.f90 \
m_lbfgs.f90 \
m_atoms.f90 \
m_ecp.f90 \
boys_function.f90 \
m_basis_set.f90 \
m_timing.f90 \
m_memory.f90 \
m_scalapack.f90 \
m_inputparam.f90 \
m_libint_tools.f90 \
m_eri.f90 m_eri_ao_mo.f90      m_eri_calculate.f90 \
m_dft_grid.f90 \
m_hamiltonian_buffer.f90 \
m_hamiltonian_onebody.f90 \
m_hamiltonian_sca.f90 \
m_hamiltonian.f90 \
io.f90



CXX_SRC = \
libint_onebody.cc  libint_onebody_grad.cc  libint_twobody.cc  libint_twobody_grad.cc

OBJ = $(SRC:.f90=.o) $(SRC:.f=.o) $(CXX_SRC:.cc=.o)

#
# Suffix rule for Fortran 90
#
%.mod :
	@if [! -f $@ ]; then \
		rm $(*F).o; \
		fi
	$(MAKE) $<

%.o : %.f90
	$(F90) $(F90_OPTS) -c -o $(*F).o $<

#
# Fortran 77 sources
# supress warning
.SUFFIXES: .o .f
.f.o:
	$(F77) -c $(F77_OPTS) $<

#
# C++ source
#
.SUFFIXES: .o .cc
.cc.o:
	$(CXX) -c $(CXX_OPTS) $<

# Libraries
lib: $(OBJ)
	ar rcs libmain.a *.o

clean:
	rm -fv libmain.a *.o *.mod


