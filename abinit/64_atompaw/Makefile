include ../make.inc

LIBNAME = lib64_atompaw.a

SRC = \
interfaces_64_atompaw.F90  m_atom.F90  m_lmn_indices.F90

include ../make.inc.post


interfaces_64_atompaw.o : interfaces_64_atompaw.F90 
m_atom.o : m_atom.F90 m_lmn_indices.o 
m_lmn_indices.o : m_lmn_indices.F90 
