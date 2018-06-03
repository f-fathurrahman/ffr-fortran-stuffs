include ../make.inc

LIBNAME = lib69_wfdesc.a

SRC = \
calc_optical_mels.F90  interfaces_69_wfdesc.F90  m_oscillators.F90  m_wfd.F90   wfd_mkrho.F90     wfd_vnlpsi.F90 \
classify_bands.F90     m_io_kss.F90              m_shirley.F90      outkss.F90  wfd_pawrhoij.F90


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



calc_optical_mels.o : calc_optical_mels.F90 m_wfd.o 
classify_bands.o : classify_bands.F90 interfaces_69_wfdesc.o m_wfd.o 
interfaces_69_wfdesc.o : interfaces_69_wfdesc.F90 m_wfd.o 
m_io_kss.o : m_io_kss.F90 
m_oscillators.o : m_oscillators.F90 
m_shirley.o : m_shirley.F90 interfaces_69_wfdesc.o m_wfd.o 
m_wfd.o : m_wfd.F90 
outkss.o : outkss.F90 m_io_kss.o 
wfd_mkrho.o : wfd_mkrho.F90 m_wfd.o 
wfd_pawrhoij.o : wfd_pawrhoij.F90 m_wfd.o 
wfd_vnlpsi.o : wfd_vnlpsi.F90 m_wfd.o 
