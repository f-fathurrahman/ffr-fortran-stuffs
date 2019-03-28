include ../make.inc

LIBNAME = lib11_qespresso_ext.a

SRC = \
flib_pwscf.F90  functionals_pwscf.F90  interfaces_11_qespresso_ext.F90 \
m_pseudo_pwscf.F90  read_upf_pwscf.F90

include ../make.inc.post

flib_pwscf.o : flib_pwscf.F90 
functionals_pwscf.o : functionals_pwscf.F90 flib_pwscf.o 
interfaces_11_qespresso_ext.o : interfaces_11_qespresso_ext.F90 
m_pseudo_pwscf.o : m_pseudo_pwscf.F90 
read_upf_pwscf.o : read_upf_pwscf.F90 interfaces_11_qespresso_ext.o m_pseudo_pwscf.o flib_pwscf.o


