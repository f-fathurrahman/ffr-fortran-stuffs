include ../make.inc

LIBNAME = lib62_poisson.a

SRC = \
interfaces_62_poisson.F90  psolver_hartree.F90  psolver_kernel.F90  psolver_rhohxc.F90


include ../make.inc.post


interfaces_62_poisson.o : interfaces_62_poisson.F90 
psolver_hartree.o : psolver_hartree.F90 interfaces_62_poisson.o 
psolver_kernel.o : psolver_kernel.F90 
psolver_rhohxc.o : psolver_rhohxc.F90 interfaces_62_poisson.o 
