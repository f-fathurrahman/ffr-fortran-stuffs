include ../make.inc

LIBNAME = lib61_occeig.a

SRC = \
dos_hdr_write.F90  init_occ_ent.F90          m_ebands.F90  newocc.F90  pareigocc.F90 \
getnel.F90         interfaces_61_occeig.F90  m_fstab.F90   occeig.F90


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



dos_hdr_write.o : dos_hdr_write.F90 
getnel.o : getnel.F90 interfaces_61_occeig.o 
init_occ_ent.o : init_occ_ent.F90 
interfaces_61_occeig.o : interfaces_61_occeig.F90 
m_ebands.o : m_ebands.F90 interfaces_61_occeig.o 
m_fstab.o : m_fstab.F90 m_ebands.o 
newocc.o : newocc.F90 interfaces_61_occeig.o 
occeig.o : occeig.F90 
pareigocc.o : pareigocc.F90 
