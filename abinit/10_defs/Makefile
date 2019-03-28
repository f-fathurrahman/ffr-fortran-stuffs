include ../make.inc

LIBNAME = lib10_defs.a

SRC = \
defs_basis.F90 \
defs_aimfields.F90  defs_fftdata.F90  defs_suscep.F90  defs_wannier90.F90 \
defs_aimprom.F90  defs_parameters.F90 defs_datatypes.F90 defs_time.F90

include ../make.inc.post

defs_aimfields.o : defs_aimfields.F90 defs_basis.o 
defs_aimprom.o : defs_aimprom.F90 defs_basis.o 
defs_basis.o : defs_basis.F90 
defs_datatypes.o : defs_datatypes.F90 defs_parameters.o defs_basis.o 
defs_fftdata.o : defs_fftdata.F90 defs_basis.o 
defs_parameters.o : defs_parameters.F90 defs_basis.o 
defs_suscep.o : defs_suscep.F90 defs_basis.o 
defs_time.o : defs_time.F90 defs_basis.o 
defs_wannier90.o : defs_wannier90.F90 defs_basis.o 
