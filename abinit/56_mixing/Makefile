include ../make.inc

LIBNAME = lib56_mixing.a

SRC = \
aprxdr.F90      dotprodm_vn.F90  interfaces_56_mixing.F90  scfcge.F90  scfopt.F90 \
dotprodm_v.F90  findminscf.F90   m_ab7_mixing.F90          scfeig.F90  sqnormm_v.F90


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



aprxdr.o : aprxdr.F90 interfaces_56_mixing.o 
dotprodm_v.o : dotprodm_v.F90 
dotprodm_vn.o : dotprodm_vn.F90 
findminscf.o : findminscf.F90 
interfaces_56_mixing.o : interfaces_56_mixing.F90 
m_ab7_mixing.o : m_ab7_mixing.F90 interfaces_56_mixing.o 
scfcge.o : scfcge.F90 interfaces_56_mixing.o 
scfeig.o : scfeig.F90 
scfopt.o : scfopt.F90 interfaces_56_mixing.o 
sqnormm_v.o : sqnormm_v.F90 
