include ../make.inc

LIBNAME = lib64_psp.a

SRC = \
atm2fft.F90         interfaces_64_psp.F90  \
m_psps.F90 \
psp1cc.F90  psp3in.F90      psp6cc.F90  psp9in.F90   vhtnzc.F90 \
cc_derivatives.F90             psp1in.F90  psp3nl.F90      psp6in.F90  pspatm.F90 \
der_int.F90         psden.F90              psp1lo.F90  psp4cc.F90      psp8cc.F90  pspcor.F90 \
dfpt_atm2fft.F90    psp10in.F90            psp1nl.F90  psp5in.F90      psp8in.F90  psp_dump_outputs.F90 \
gg1cc.F90           psp10nl.F90            psp2in.F90  psp5lo.F90      psp8lo.F90  pspini.F90 \
gp1cc.F90           psp11lo.F90            psp2lo.F90  psp5nl.F90      psp8nl.F90  sincos.F90 \
gpp1cc.F90          psp11nl.F90            psp2nl.F90  psp6cc_drh.F90  upf2abinit.F90

# not included: psp9cc.F90

include ../make.inc.post


atm2fft.o : atm2fft.F90 
cc_derivatives.o : cc_derivatives.F90 
der_int.o : der_int.F90 
dfpt_atm2fft.o : dfpt_atm2fft.F90 
gg1cc.o : gg1cc.F90 
gp1cc.o : gp1cc.F90 
gpp1cc.o : gpp1cc.F90 
interfaces_64_psp.o : interfaces_64_psp.F90 
m_psps.o : m_psps.F90 
psden.o : psden.F90 
psp10in.o : psp10in.F90 interfaces_64_psp.o 
psp10nl.o : psp10nl.F90 
psp11lo.o : psp11lo.F90 
psp11nl.o : psp11nl.F90 
psp1cc.o : psp1cc.F90 interfaces_64_psp.o 
psp1in.o : psp1in.F90 interfaces_64_psp.o 
psp1lo.o : psp1lo.F90 interfaces_64_psp.o 
psp1nl.o : psp1nl.F90 interfaces_64_psp.o 
psp2in.o : psp2in.F90 interfaces_64_psp.o 
psp2lo.o : psp2lo.F90 
psp2nl.o : psp2nl.F90 
psp3in.o : psp3in.F90 interfaces_64_psp.o 
psp3nl.o : psp3nl.F90 
psp4cc.o : psp4cc.F90 
psp5in.o : psp5in.F90 interfaces_64_psp.o 
psp5lo.o : psp5lo.F90 
psp5nl.o : psp5nl.F90 
psp6cc_drh.o : psp6cc_drh.F90 interfaces_64_psp.o 
psp6cc.o : psp6cc.F90 interfaces_64_psp.o 
psp6in.o : psp6in.F90 interfaces_64_psp.o 
psp8cc.o : psp8cc.F90 
psp8in.o : psp8in.F90 interfaces_64_psp.o m_psps.o 
psp8lo.o : psp8lo.F90 
psp8nl.o : psp8nl.F90 
#psp9cc.o : psp9cc.F90 
psp9in.o : psp9in.F90 interfaces_64_psp.o m_psps.o 
pspatm.o : pspatm.F90 interfaces_64_psp.o m_psps.o 
pspcor.o : pspcor.F90 
psp_dump_outputs.o : psp_dump_outputs.F90 
pspini.o : pspini.F90 interfaces_64_psp.o m_psps.o 
sincos.o : sincos.F90 
upf2abinit.o : upf2abinit.F90 interfaces_64_psp.o 
vhtnzc.o : vhtnzc.F90 
