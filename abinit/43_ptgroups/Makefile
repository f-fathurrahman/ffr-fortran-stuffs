include ../make.inc

LIBNAME = lib43_ptgroups.a

SRC = \
interfaces_43_ptgroups.F90  ptg_C2h.F90  ptg_C3v.F90  ptg_C6h.F90  ptg_D2.F90   ptg_D4.F90   ptg_Oh.F90 \
m_defs_ptgroups.F90         ptg_C2v.F90  ptg_C4.F90   ptg_C6v.F90  ptg_D2h.F90  ptg_D4h.F90  ptg_S4.F90 \
m_ptgroups.F90              ptg_C3.F90   ptg_C4h.F90  ptg_Ci.F90   ptg_D3d.F90  ptg_D6.F90   ptg_Td.F90 \
ptg_C1.F90                  ptg_C3h.F90  ptg_C4v.F90  ptg_Cs.F90   ptg_D3.F90   ptg_D6h.F90  ptg_T.F90 \
ptg_C2.F90                  ptg_C3i.F90  ptg_C6.F90   ptg_D2d.F90  ptg_D3h.F90  ptg_O.F90    ptg_Th.F90


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



interfaces_43_ptgroups.o : interfaces_43_ptgroups.F90 m_defs_ptgroups.o 
m_defs_ptgroups.o : m_defs_ptgroups.F90 
m_ptgroups.o : m_ptgroups.F90 interfaces_43_ptgroups.o m_defs_ptgroups.o 
ptg_C1.o : ptg_C1.F90 m_defs_ptgroups.o 
ptg_C2.o : ptg_C2.F90 m_defs_ptgroups.o 
ptg_C2h.o : ptg_C2h.F90 m_defs_ptgroups.o 
ptg_C2v.o : ptg_C2v.F90 m_defs_ptgroups.o 
ptg_C3.o : ptg_C3.F90 m_defs_ptgroups.o 
ptg_C3h.o : ptg_C3h.F90 m_defs_ptgroups.o 
ptg_C3i.o : ptg_C3i.F90 m_defs_ptgroups.o 
ptg_C3v.o : ptg_C3v.F90 m_defs_ptgroups.o 
ptg_C4.o : ptg_C4.F90 m_defs_ptgroups.o 
ptg_C4h.o : ptg_C4h.F90 m_defs_ptgroups.o 
ptg_C4v.o : ptg_C4v.F90 m_defs_ptgroups.o 
ptg_C6.o : ptg_C6.F90 m_defs_ptgroups.o 
ptg_C6h.o : ptg_C6h.F90 m_defs_ptgroups.o 
ptg_C6v.o : ptg_C6v.F90 m_defs_ptgroups.o 
ptg_Ci.o : ptg_Ci.F90 m_defs_ptgroups.o 
ptg_Cs.o : ptg_Cs.F90 m_defs_ptgroups.o 
ptg_D2d.o : ptg_D2d.F90 m_defs_ptgroups.o 
ptg_D2.o : ptg_D2.F90 m_defs_ptgroups.o 
ptg_D2h.o : ptg_D2h.F90 m_defs_ptgroups.o 
ptg_D3d.o : ptg_D3d.F90 m_defs_ptgroups.o 
ptg_D3.o : ptg_D3.F90 m_defs_ptgroups.o 
ptg_D3h.o : ptg_D3h.F90 m_defs_ptgroups.o 
ptg_D4.o : ptg_D4.F90 m_defs_ptgroups.o 
ptg_D4h.o : ptg_D4h.F90 m_defs_ptgroups.o 
ptg_D6.o : ptg_D6.F90 m_defs_ptgroups.o 
ptg_D6h.o : ptg_D6h.F90 m_defs_ptgroups.o 
ptg_O.o : ptg_O.F90 m_defs_ptgroups.o 
ptg_Oh.o : ptg_Oh.F90 m_defs_ptgroups.o 
ptg_S4.o : ptg_S4.F90 m_defs_ptgroups.o 
ptg_Td.o : ptg_Td.F90 m_defs_ptgroups.o 
ptg_T.o : ptg_T.F90 m_defs_ptgroups.o 
ptg_Th.o : ptg_Th.F90 m_defs_ptgroups.o 
