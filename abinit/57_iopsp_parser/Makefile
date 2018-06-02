include ../make.inc

LIBNAME = lib57_iopsp_parser.a

SRC = \
m_psxml2ab.F90 \
inpspheads.F90  \
interfaces_57_iopsp_parser.F90  pawpsxml2ab.F90  upfheader2abi.F90

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


inpspheads.o : inpspheads.F90 interfaces_57_iopsp_parser.o m_psxml2ab.o 
interfaces_57_iopsp_parser.o : interfaces_57_iopsp_parser.F90 
m_psxml2ab.o : m_psxml2ab.F90 
pawpsxml2ab.o : pawpsxml2ab.F90 
upfheader2abi.o : upfheader2abi.F90 interfaces_57_iopsp_parser.o

