include ../make.inc

LIBNAME = lib59_ionetcdf.a

SRC = \
interfaces_59_ionetcdf.F90  m_crystal_io.F90  m_dfpt_io.F90 \
write_eig.F90  wrt_moldyn_netcdf.F90


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



interfaces_59_ionetcdf.o : interfaces_59_ionetcdf.F90 
m_crystal_io.o : m_crystal_io.F90 
m_dfpt_io.o : m_dfpt_io.F90 
write_eig.o : write_eig.F90 
wrt_moldyn_netcdf.o : wrt_moldyn_netcdf.F90 
