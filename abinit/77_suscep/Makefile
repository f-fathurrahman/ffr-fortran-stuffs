include ../make.inc

LIBNAME = lib77_suscep.a

SRC = \
interfaces_77_suscep.F90  suscep_stat.F90  susk.F90  suskmm.F90


include ../make.inc.post


interfaces_77_suscep.o : interfaces_77_suscep.F90 
suscep_stat.o : suscep_stat.F90 interfaces_77_suscep.o 
susk.o : susk.F90 
suskmm.o : suskmm.F90 
