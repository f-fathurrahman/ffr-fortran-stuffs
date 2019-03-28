include ../make.inc

LIBNAME = lib68_rsprc.a

SRC = \
interfaces_68_rsprc.F90  mod_prc_memory.F90  newvtr.F90  prcref_PMA.F90    prcrskerker2.F90 \
moddiel_csrb.F90         newrho.F90          prcref.F90  prcrskerker1.F90  wvl_prcref.F90


include ../make.inc.post


interfaces_68_rsprc.o : interfaces_68_rsprc.F90 
moddiel_csrb.o : moddiel_csrb.F90 
mod_prc_memory.o : mod_prc_memory.F90 
newrho.o : newrho.F90 interfaces_68_rsprc.o 
newvtr.o : newvtr.F90 interfaces_68_rsprc.o 
prcref.o : prcref.F90 interfaces_68_rsprc.o 
prcref_PMA.o : prcref_PMA.F90 interfaces_68_rsprc.o 
prcrskerker1.o : prcrskerker1.F90 mod_prc_memory.o 
prcrskerker2.o : prcrskerker2.F90 
wvl_prcref.o : wvl_prcref.F90 
