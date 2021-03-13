F90 = gfortran
F90FLAGS = -Wall -fcheck=bounds #-O3

PROG = libmain.a

OBJS = stop_error.o \
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
schroed_inward_adams.o 

all: $(PROG)

$(PROG): $(OBJS)
	ar rcs $@ $(OBJS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

