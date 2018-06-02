include ../make.inc

LIBNAME = lib02_clib.a

SRC = \
interfaces_02_clib.F90  m_clib.F90

C_SRC = \
align.c           cclock.c  fsi_posix.c  intrinsics.c  md5.c           xexit.c \
calling_levmar.c  etime.c   gnu_tools.c  mallinfo.c    progress_bar.c  xmalloc.c


OBJ = $(SRC:.F90=.o) $(C_SRC:.c=.o)

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



