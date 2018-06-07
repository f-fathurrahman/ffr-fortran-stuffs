include make.inc

LIBS_EXT =


SRC = \
ank.f90                                    m_mpi_world.f90 \
basis_function_dipole.f90                  gos_basis_function.f90        m_scalapack.f90 \
basis_function_prod.f90                    header.f90                    m_timing.f90 \
basis_function_quadrupole.f90              init_basis_function.f90       m_tools.f90 \
boys_function.f90                          init_basis_function_pure.f90  mulliken_pdos.f90 \
calculate_basis_functions_gradr_batch.f90  init_basis_set.f90            m_warning.f90 \
calculate_basis_functions_gradr.f90        init_dft_grid.f90             nuclear_basis_function.f90 \
calculate_basis_functions_laplr.f90        init_gaussian_general.f90     nucleus_recurrence.f90 \
calculate_basis_functions_r_batch.f90      kinetic_basis_function.f90    number_basis_function_am.f90 \
calculate_basis_functions_r.f90            kinetic_recurrence.f90        numerical_kinetic.f90 \
cnk.f90                                    m_atoms.f90                   numerical_nucleus.f90 \
compare_basis_function.f90                 m_basis_set.f90               output_new_homolumo.f90 \
compare_basis_set.f90                      m_cart_to_pure.f90            overlap_basis_function.f90 \
compare_gaussian.f90                       m_definitions.f90             overlap_recurrence.f90 \
destroy_basis_function.f90                 m_dft_grid.f90                overlap_three_basis_function.f90 \
destroy_basis_set.f90                      m_ecp.f90                     plot_cube_wfn.f90 \
destroy_cart_to_pure_transforms.f90        m_elements.f90                plot_rho.f90 \
destroy_dft_grid.f90                       m_eri_ao_mo.f90               plot_rho_list.f90 \
double_factorial.f90                       m_eri_calculate.f90           plot_wfn.f90 \
dump_out_energy.f90                        m_eri.f90                     prepare_basis_functions_gradr.f90 \
dump_out_matrix.f90                        m_gaussian.f90                prepare_basis_functions_r.f90 \
dump_out_occupation.f90                    m_gos.f90                     print_basis_function.f90 \
echo_input_variables.f90                   m_hamiltonian_buffer.f90      print_gaussian.f90 \
eval_basis_function.f90                    m_hamiltonian.f90             product_gaussian.f90 \
eval_basis_function_grad.f90               m_hamiltonian_onebody.f90     read_basis_function.f90 \
eval_basis_function_lapl.f90               m_hamiltonian_sca.f90         read_basis_set.f90 \
eval_gaussian.f90                          m_inputparam.f90              read_energy_qp.f90 \
eval_gaussian_grad.f90                     m_lbfgs.f90                   setup_bf_radius.f90 \
eval_gaussian_lapl.f90                     m_libint_tools.f90            setup_cart_to_pure_transforms.f90 \
evaluate_gos.f90                           m_memory.f90                  smooth_step.f90 \
evaluate_wfn_r.f90                         m_mpi_auxil.f90               wfn_parity.f90 \
get_basis_functions_gradr_batch.f90        m_mpi.f90                     wfn_reflection.f90 \
get_basis_functions_gradr.f90              m_mpi_grid.f90                write_basis_function.f90 \
get_basis_functions_r_batch.f90            m_mpi_local.f90               write_basis_set.f90 \
get_basis_functions_r.f90                  m_mpi_ortho.f90               write_energy_qp.f90 \
get_gaussian_type_tag.f90                  m_mpi_trans.f90




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


include make.depends
