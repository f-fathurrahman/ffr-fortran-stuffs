include ../make.inc

LIBNAME = lib10_defs.a

SRC = \
defs_rectypes.F90  first_rec.F90     interfaces_68_recursion.F90  m_rec_tools.F90       recursion.F90 \
density_rec.F90    getngrec.F90      m_hidecudarec.F90            nlenergyrec.F90       recursion_nl.F90 \
entropyrec.F90     gran_potrec.F90   m_pretty_rec.F90             pspnl_hgh_rec.F90     vn_nl_rec.F90 \
fermisolverec.F90  green_kernel.F90  m_rec.F90                    pspnl_operat_rec.F90  vtorhorec.F90



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



defs_rectypes.o : defs_rectypes.F90 
density_rec.o : density_rec.F90 m_rec_tools.o 
entropyrec.o : entropyrec.F90 
fermisolverec.o : fermisolverec.F90 interfaces_68_recursion.o defs_rectypes.o 
first_rec.o : first_rec.F90 interfaces_68_recursion.o m_hidecudarec.o m_rec_tools.o m_rec.o defs_rectypes.o 
getngrec.o : getngrec.F90 
gran_potrec.o : gran_potrec.F90 
green_kernel.o : green_kernel.F90 
interfaces_68_recursion.o : interfaces_68_recursion.F90 defs_rectypes.o 
m_hidecudarec.o : m_hidecudarec.F90 defs_rectypes.o 
m_pretty_rec.o : m_pretty_rec.F90 
m_rec.o : m_rec.F90 interfaces_68_recursion.o m_rec_tools.o defs_rectypes.o 
m_rec_tools.o : m_rec_tools.F90 defs_rectypes.o 
nlenergyrec.o : nlenergyrec.F90 interfaces_68_recursion.o m_rec_tools.o defs_rectypes.o 
pspnl_hgh_rec.o : pspnl_hgh_rec.F90 defs_rectypes.o 
pspnl_operat_rec.o : pspnl_operat_rec.F90 defs_rectypes.o 
recursion.o : recursion.F90 interfaces_68_recursion.o m_rec_tools.o defs_rectypes.o 
recursion_nl.o : recursion_nl.F90 interfaces_68_recursion.o m_rec_tools.o defs_rectypes.o 
vn_nl_rec.o : vn_nl_rec.F90 defs_rectypes.o 
vtorhorec.o : vtorhorec.F90 interfaces_68_recursion.o m_rec_tools.o m_rec.o m_pretty_rec.o defs_rectypes.o 
