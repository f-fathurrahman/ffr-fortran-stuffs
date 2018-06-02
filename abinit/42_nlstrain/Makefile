include ../make.inc

LIBNAME = lib42_nlstrain.a

SRC = \
contistr01.F90  contistr12.F90  contstr23.F90   contstr25.F90  interfaces_42_nlstrain.F90 \
contistr03.F90  contstr21.F90   contstr25a.F90  contstr26.F90


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



contistr01.o : contistr01.F90 
contistr03.o : contistr03.F90 
contistr12.o : contistr12.F90 
contstr21.o : contstr21.F90 
contstr23.o : contstr23.F90 
contstr25a.o : contstr25a.F90 
contstr25.o : contstr25.F90 
contstr26.o : contstr26.F90 
interfaces_42_nlstrain.o : interfaces_42_nlstrain.F90 
