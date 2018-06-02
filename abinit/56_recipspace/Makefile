include ../make.inc

LIBNAME = lib56_recipspace.a

SRC = \
interfaces_56_recipspace.F90 \
cg_rotate.F90       get_kpt_fullbz.F90            irrzg.F90      m_ab7_kpoints.F90  m_kpts.F90     setsym.F90 \
getcut.F90          getmpw.F90                    kpgio.F90      m_bz_mesh.F90      m_nesting.F90  smpbz.F90 \
get_full_kgrid.F90  getph.F90                     kpgstr.F90     m_double_grid.F90 \
m_gsphere.F90 \
m_skw.F90      symg.F90 \
getkgrid.F90        initylmg.F90                  laplacian.F90  ph1d3d.F90     symkchk.F90 \
getkpgnorm.F90  listkk.F90     mkkin.F90          setmqgrid.F90  testkgrid.F90


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



