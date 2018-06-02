include ../make.inc

LIBNAME = lib57_iovars.a

SRC = \
interfaces_57_iovars.F90 \
append_xyz.F90  compute_kgb_indicator.F90  inqpt.F90                 macroin.F90      outvar_i_n.F90 \
chkdpr.F90      create_nc_file.F90         memana.F90       outvar_o_z.F90 \
chkinp.F90      finddistrproc.F90          invacuum.F90              memorf.F90       outvars.F90 \
chkint_eq.F90   getdim_nloc.F90            invars0.F90               memory_eval.F90  outxml.F90 \
chkint.F90      importxyz.F90              invars1.F90               memory.F90       parsefile.F90 \
chkint_ge.F90   indefo1.F90                invars1m.F90              m_ingeo_img.F90  prtocc.F90 \
chkint_le.F90   indefo.F90                 invars2.F90               mpi_setup.F90    prttagm.F90 \
chkint_ne.F90   ingeobld.F90               invars2m.F90              out_acknowl.F90  prttagm_images.F90 \
chkint_prt.F90  ingeo.F90                  m_ab7_invars_f90.F90      outqmc.F90       write_var_netcdf.F90 \
chkvars.F90     inkpts.F90                 macroin2.F90              outvar_a_h.F90   wvl_memory.F90


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



