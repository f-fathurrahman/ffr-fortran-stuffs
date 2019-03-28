include ../make.inc

LIBNAME = lib94_scfcv.a

SRC = \
afterscfloop.F90  elpolariz.F90  interfaces_94_scfcv.F90  m_scfcv.F90  outscfcv.F90  scfcv.F90


include ../make.inc.post


afterscfloop.o : afterscfloop.F90 interfaces_94_scfcv.o 
elpolariz.o : elpolariz.F90 
interfaces_94_scfcv.o : interfaces_94_scfcv.F90 
m_scfcv.o : m_scfcv.F90 interfaces_94_scfcv.o 
outscfcv.o : outscfcv.F90 
scfcv.o : scfcv.F90 interfaces_94_scfcv.o 
