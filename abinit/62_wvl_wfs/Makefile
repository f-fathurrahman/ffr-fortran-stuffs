include ../make.inc

LIBNAME = lib62_wvl_defs.a

SRC = \
interfaces_62_wvl_wfs.F90  wvl_nl_gradient.F90  wvl_rwwf.F90 \
wvl_hpsitopsi.F90          wvl_psitohpsi.F90    wvl_tail_corrections.F90


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



interfaces_62_wvl_wfs.o : interfaces_62_wvl_wfs.F90 
wvl_hpsitopsi.o : wvl_hpsitopsi.F90 
wvl_nl_gradient.o : wvl_nl_gradient.F90 
wvl_psitohpsi.o : wvl_psitohpsi.F90 
wvl_rwwf.o : wvl_rwwf.F90 
wvl_tail_corrections.o : wvl_tail_corrections.F90 
