include ../make.inc

LIBNAME = lib62_ctqmc.a

SRC = \
interfaces_62_ctqmc.F90  m_CtqmcInterface.F90    m_global.F90 \
m_ListCdagC.F90  m_OurRng.F90  m_VectorInt.F90 \
m_BathOperator.F90       m_data4entropyDMFT.F90  m_GreenHyb.F90          m_MapHyb.F90     m_Stat.F90 \
m_Ctqmc.F90              m_FFTHyb.F90            m_ImpurityOperator.F90  m_MatrixHyb.F90  m_Vector.F90


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



interfaces_62_ctqmc.o : interfaces_62_ctqmc.F90 
m_BathOperator.o : m_BathOperator.F90 BathOperator_hybrid BathOperator_hybrid.h defs.h m_ListCdagC.o m_global.o m_VectorInt.o m_Vector.o m_MatrixHyb.o 
m_Ctqmc.o : m_Ctqmc.F90 defs.h m_Vector.o m_OurRng.o m_FFTHyb.o m_Stat.o m_ImpurityOperator.o m_BathOperator.o m_GreenHyb.o m_global.o 
m_CtqmcInterface.o : m_CtqmcInterface.F90 defs.h m_Ctqmc.o 
m_data4entropyDMFT.o : m_data4entropyDMFT.F90 
m_FFTHyb.o : m_FFTHyb.F90 defs.h m_global.o 
m_global.o : m_global.F90 defs.h 
m_GreenHyb.o : m_GreenHyb.F90 defs.h m_MapHyb.o m_ListCdagC.o m_VectorInt.o m_Vector.o m_MatrixHyb.o m_global.o 
m_ImpurityOperator.o : m_ImpurityOperator.F90 ListCdagC_firstHigher ListCdagC_firstHigher.h defs.h m_global.o m_ListCdagC.o 
m_ListCdagC.o : m_ListCdagC.F90 ListCdagC_firstHigher ListCdagC_firstHigher.h defs.h m_global.o 
m_MapHyb.o : m_MapHyb.F90 defs.h m_global.o 
m_MatrixHyb.o : m_MatrixHyb.F90 defs.h m_global.o 
m_OurRng.o : m_OurRng.F90 defs.h 
m_Stat.o : m_Stat.F90 defs.h m_global.o 
m_Vector.o : m_Vector.F90 defs.h m_global.o 
m_VectorInt.o : m_VectorInt.F90 defs.h m_global.o 
