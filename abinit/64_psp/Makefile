include ../make.inc

LIBNAME = lib64_psp.a

SRC = \
atm2fft.F90         interfaces_64_psp.F90  psp1cc.F90  psp3in.F90      psp6cc.F90  psp9in.F90   vhtnzc.F90 \
cc_derivatives.F90  m_psps.F90             psp1in.F90  psp3nl.F90      psp6in.F90  pspatm.F90 \
der_int.F90         psden.F90              psp1lo.F90  psp4cc.F90      psp8cc.F90  pspcor.F90 \
dfpt_atm2fft.F90    psp10in.F90            psp1nl.F90  psp5in.F90      psp8in.F90  psp_dump_outputs.F90 \
gg1cc.F90           psp10nl.F90            psp2in.F90  psp5lo.F90      psp8lo.F90  pspini.F90 \
gp1cc.F90           psp11lo.F90            psp2lo.F90  psp5nl.F90      psp8nl.F90  sincos.F90 \
gpp1cc.F90          psp11nl.F90            psp2nl.F90  psp6cc_drh.F90  psp9cc.F90  upf2abinit.F90


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



