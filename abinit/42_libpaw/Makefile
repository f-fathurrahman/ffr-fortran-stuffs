include ../make.inc

LIBNAME = lib42_libpaw.a

SRC = \
interfaces_42_libpaw.F90  m_libpaw_tools.F90 \
m_paw_numeric.F90 \
m_pawrad.F90 \
m_pawtab.F90 \
m_paw_atom.F90 \
m_paw_sphharm.F90 \
m_paw_finegrid.F90 \
m_libpaw_defs.F90  m_paral_atom.F90    m_pawcprj.F90    m_paw_gaussfit.F90 \
m_libpaw_libxc.F90 \
m_pawang.F90 \
m_pawxmlps.F90 \
m_pawxc.F90 \
m_pawpsp.F90 \
m_paw_an.F90 \
m_paw_io.F90 \
m_paw_ij.F90 \
m_pawfgrtab.F90 \
m_pawrhoij.F90 \
m_pawdij.F90 \
m_libpaw_mpi.F90        


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



