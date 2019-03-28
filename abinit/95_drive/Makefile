include ../make.inc

LIBNAME = lib95_drive.a

SRC = \
bethe_salpeter.F90  dtfil_init.F90       gstate.F90               iofn1.F90         nonlinear.F90    sigma.F90 \
dfpt_looppert.F90   dtfil_init_img.F90   gstateimg.F90            m_argparse.F90    pawuj_drive.F90  testfi.F90 \
dfptnl_loop.F90     dtfil_init_time.F90  gwls_sternheimer.F90     mover_effpot.F90  respfn.F90       timana.F90 \
driver.F90          eph.F90              interfaces_95_drive.F90  mover.F90         screening.F90    wfk_analyze.F90

include ../make.inc.post


bethe_salpeter.o : bethe_salpeter.F90 
dfpt_looppert.o : dfpt_looppert.F90 
dfptnl_loop.o : dfptnl_loop.F90 
driver.o : driver.F90 interfaces_95_drive.o 
dtfil_init.o : dtfil_init.F90 
dtfil_init_img.o : dtfil_init_img.F90 
dtfil_init_time.o : dtfil_init_time.F90 
eph.o : eph.F90 
gstate.o : gstate.F90 interfaces_95_drive.o 
gstateimg.o : gstateimg.F90 interfaces_95_drive.o 
gwls_sternheimer.o : gwls_sternheimer.F90 interfaces_95_drive.o 
interfaces_95_drive.o : interfaces_95_drive.F90 
iofn1.o : iofn1.F90 
m_argparse.o : m_argparse.F90 
mover_effpot.o : mover_effpot.F90 interfaces_95_drive.o 
mover.o : mover.F90 interfaces_95_drive.o 
nonlinear.o : nonlinear.F90 interfaces_95_drive.o 
pawuj_drive.o : pawuj_drive.F90 
respfn.o : respfn.F90 interfaces_95_drive.o 
screening.o : screening.F90 
sigma.o : sigma.F90 
testfi.o : testfi.F90 
timana.o : timana.F90 
wfk_analyze.o : wfk_analyze.F90 


