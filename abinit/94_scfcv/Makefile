include ../make.inc

LIBNAME = lib94_scfcv.a

SRC = \
afterscfloop.F90  elpolariz.F90  interfaces_94_scfcv.F90  m_scfcv.F90  outscfcv.F90  scfcv.F90


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



afterscfloop.o : afterscfloop.F90 interfaces_94_scfcv.o 
elpolariz.o : elpolariz.F90 
interfaces_94_scfcv.o : interfaces_94_scfcv.F90 
m_scfcv.o : m_scfcv.F90 interfaces_94_scfcv.o 
outscfcv.o : outscfcv.F90 
scfcv.o : scfcv.F90 interfaces_94_scfcv.o 
