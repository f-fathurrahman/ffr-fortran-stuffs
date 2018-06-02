include ../make.inc

LIBNAME = lib45_geomoptim.a

SRC = \
m_abimover.F90 \
m_pimd.F90 \
m_abihist.F90 \
m_bfgs.F90 \
m_mep.F90 \
m_ga.F90 \
interfaces_45_geomoptim.F90 \
calc_b_matrix.F90   pred_bfgs.F90  pred_isokinetic.F90  pred_verlet.F90 \
deloc2xcart.F90   m_lbfgs.F90  pred_delocint.F90  pred_isothermal.F90  prtxfase.F90 \
fred2fdeloc.F90   pred_diisrelax.F90    pred_langevin.F90    strainsym.F90 \
pred_hmc.F90    pred_lbfgs.F90   xcart2deloc.F90 \
isotemp.F90    m_use_ga.F90   predict_copy.F90      pred_moldyn.F90      xfh_recover_deloc.F90 \
pimd_langevin_npt.F90    predictimg.F90        pred_nose.F90        xfh_recover_new.F90 \
m_abiimages.F90  pimd_langevin_nvt.F90    predict_neb.F90       pred_simple.F90      xfh_update.F90 \
pimd_nosehoover_npt.F90  predict_pimd.F90      pred_srkna14.F90     xfpack_f2vout.F90 \
matpointsym.F90  pimd_nosehoover_nvt.F90  predict_steepest.F90  pred_steepdesc.F90   xfpack_vin2x.F90 \
prec_simple.F90  predict_string.F90    pred_velverlet.F90   xfpack_x2vin.F90


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



