include ../make.inc

LIBNAME = lib52_ffr_mpi_noabirule.a

SRC = \
interfaces_52_fft_mpi_noabirule.F90 \
m_fftcore.F90  \
m_fft_mesh.F90 \
m_fftw3.F90 m_sgfft.F90 \
m_dfti.F90  m_sg2002.F90 \
bound.F90  indirect_parallel_Fourier.F90  sphereboundary.F90 \
getng.F90


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



