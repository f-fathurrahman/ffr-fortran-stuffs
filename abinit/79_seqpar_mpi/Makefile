include ../make.inc

LIBNAME = lib79_seqpar_mpi.a

SRC = \
dfptnl_mv.F90                 lobpcgwf.F90    tddft.F90   wfsinp.F90   wvl_wfsinp_scratch.F90 \
interfaces_79_seqpar_mpi.F90  m_lobpcg.F90    vtorho.F90  wvl_wfsinp_disk.F90 \
inwffil.F90                   m_lobpcgwf.F90  vtowfk.F90  wvl_wfsinp_reformat.F90


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



dfptnl_mv.o : dfptnl_mv.F90 
interfaces_79_seqpar_mpi.o : interfaces_79_seqpar_mpi.F90 
inwffil.o : inwffil.F90 interfaces_79_seqpar_mpi.o 
lobpcgwf.o : lobpcgwf.F90 m_lobpcg.o 
m_lobpcg.o : m_lobpcg.F90 
m_lobpcgwf.o : m_lobpcgwf.F90 m_lobpcg.o 
tddft.o : tddft.F90 
vtorho.o : vtorho.F90 interfaces_79_seqpar_mpi.o 
vtowfk.o : vtowfk.F90 interfaces_79_seqpar_mpi.o m_lobpcgwf.o 
wfsinp.o : wfsinp.F90 
wvl_wfsinp_disk.o : wvl_wfsinp_disk.F90 
wvl_wfsinp_reformat.o : wvl_wfsinp_reformat.F90 
wvl_wfsinp_scratch.o : wvl_wfsinp_scratch.F90 
