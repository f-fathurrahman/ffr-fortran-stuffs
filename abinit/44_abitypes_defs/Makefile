include ../make.inc

LIBNAME = lib44_abitypes.a

SRC = \
m_distribfft.F90 \
defs_abitypes.F90  m_args_gs.F90 m_efield.F90  m_energies.F90 \
m_results_gs.F90 \
m_results_img.F90 \
interfaces_44_abitypes_defs.F90  m_efmas_defs.F90  m_results_respfn.F90

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



