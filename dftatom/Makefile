F90 = gfortran
F90FLAGS = -Wall -fcheck=bounds #-O3

PROG = libmain.a

OBJS = \
m_dft_data.o \
m_states.o \
stop_error.o \
gen_rmesh_exp.o gen_drmesh_exp.o \
integrate_trapz.o \
rk4_integrate.o \
integrate_rschroed_rk4.o \
integrate_rproblem_outward.o \
adams_interp_outward_implicit.o \
schroed_outward_adams.o \
calc_interp.o \
get_interp_val.o \
get_midpoints.o \
get_min_idx.o \
get_n_nodes.o \
solve_radial_eigenproblem.o \
find_ctp.o \
is_E_above.o \
integrate_rproblem_inward.o \
adams_interp_inward_implicit.o \
schroed_inward_adams.o \
adams_interp_outward.o \
adams_extrapolation_outward.o \
rpoisson_outward_pc.o \
rpoisson_outward_rk4.o \
rk4_integrate3.o \
get_n_orb.o \
get_hydrogen_energies.o \
calc_E_nl.o \
calc_TF_potential.o \
atom_lda.o \
scf_potmix_anderson.o \
V_to_rho.o \
calc_Vxc.o \
calc_Vh.o \
calc_Etot.o \
rho_to_V.o \
ks_potential_map.o

all: $(PROG)

$(PROG): $(OBJS)
	ar rcs $@ $(OBJS)

clean:
	rm -f $(PROG) *.o *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

