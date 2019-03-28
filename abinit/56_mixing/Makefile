include ../make.inc

LIBNAME = lib56_mixing.a

SRC = \
aprxdr.F90      dotprodm_vn.F90  interfaces_56_mixing.F90  scfcge.F90  scfopt.F90 \
dotprodm_v.F90  findminscf.F90   m_ab7_mixing.F90          scfeig.F90  sqnormm_v.F90


include ../make.inc.post


aprxdr.o : aprxdr.F90 interfaces_56_mixing.o 
dotprodm_v.o : dotprodm_v.F90 
dotprodm_vn.o : dotprodm_vn.F90 
findminscf.o : findminscf.F90 
interfaces_56_mixing.o : interfaces_56_mixing.F90 
m_ab7_mixing.o : m_ab7_mixing.F90 interfaces_56_mixing.o 
scfcge.o : scfcge.F90 interfaces_56_mixing.o 
scfeig.o : scfeig.F90 
scfopt.o : scfopt.F90 interfaces_56_mixing.o 
sqnormm_v.o : sqnormm_v.F90 
