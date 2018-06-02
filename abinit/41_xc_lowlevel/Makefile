include ../make.inc

LIBNAME = lib41_xc_lowlevel.a

SRC = \
interfaces_41_xc_lowlevel.F90 \
m_libxc_functionals.F90 \
check_kxc.F90 \
size_dvxc.F90  xclb.F90        xcpzca.F90  xcwign.F90 \
drivexc.F90       invcb.F90   xchcth.F90     xcmult.F90      xcspol.F90  xcxalp.F90 \
drivexc_main.F90  mkdenpos.F90 xchelu.F90     xcpbe.F90       xctetr.F90 \
echo_xc_name.F90  xciit.F90      xcpositron.F90  xctfw.F90


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



