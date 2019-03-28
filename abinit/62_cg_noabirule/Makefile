include ../make.inc

LIBNAME = lib62_cg_noabirule.a

SRC = \
bracketing.F90  cgpr.F90        dummy_functions.inc  frskerker2.F90  linmin.F90 \
brent.F90       dotproduct.F90  frskerker1.F90       interfaces_62_cg_noabirule.F90

include ../make.inc.post

bracketing.o : bracketing.F90 dummy_functions.inc 
brent.o : brent.F90 dummy_functions.inc interfaces_62_cg_noabirule.o 
cgpr.o : cgpr.F90 dummy_functions.inc interfaces_62_cg_noabirule.o 
dotproduct.o : dotproduct.F90 
frskerker1.o : frskerker1.F90 interfaces_62_cg_noabirule.o 
frskerker2.o : frskerker2.F90 interfaces_62_cg_noabirule.o 
interfaces_62_cg_noabirule.o : interfaces_62_cg_noabirule.F90 
linmin.o : linmin.F90 dummy_functions.inc interfaces_62_cg_noabirule.o 
