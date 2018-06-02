include ../make.inc

LIBNAME = lib51_manage_mpi.a

SRC = \
interfaces_51_manage_mpi.F90 \
distrb2.F90        initmpi_atom.F90  initmpi_pert.F90              m_io_redirect.F90    pspheads_comm.F90 \
distrb2_hf.F90     initmpi_band.F90 \
m_mpinfo.F90 \
initmpi_seq.F90 \
get_npert_rbz.F90  initmpi_grid.F90  initmpi_world.F90             m_wffile.F90 \
herald.F90         initmpi_img.F90   m_xredistribute.F90


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



