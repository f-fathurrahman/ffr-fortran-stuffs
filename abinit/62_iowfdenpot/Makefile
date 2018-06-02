include ../make.inc

LIBNAME = lib62_iowfdenpot.a

SRC = \
initwf.F90                    m_io_screening.F90  out1dm.F90   randac.F90  WffReadEigK.F90 \
interfaces_62_iowfdenpot.F90  m_wfk.F90           outwant.F90  rdnpw.F90   WffReadSkipK.F90


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



initwf.o : initwf.F90 interfaces_62_iowfdenpot.o m_wfk.o 
interfaces_62_iowfdenpot.o : interfaces_62_iowfdenpot.F90 
m_io_screening.o : m_io_screening.F90 
m_wfk.o : m_wfk.F90 
out1dm.o : out1dm.F90 
outwant.o : outwant.F90 
randac.o : randac.F90 
rdnpw.o : rdnpw.F90 
WffReadEigK.o : WffReadEigK.F90 
WffReadSkipK.o : WffReadSkipK.F90 
