include ../make.inc

LIBNAME = lib66_nonlocal.a

SRC = \
cont13.F90     cont35.F90                  metcon.F90         mkffkg.F90      nonlop_test.F90  opernlb_ylm.F90 \
cont22cso.F90  cont3.F90                   metcon_so.F90      mkffnl.F90      nonlop_ylm.F90   opernlc_ylm.F90 \
cont22.F90     ctocprj.F90                 metric_so.F90      mkkpg.F90       opernl2.F90      opernld_ylm.F90 \
cont22so.F90   ddkten.F90                  metstr.F90         m_vkbr.F90      opernl3.F90      scalewf_nonlop.F90 \
cont24.F90     dfpt_mkffkg.F90             m_fock.F90         nonlop.F90      opernl4a.F90     strsocv.F90 \
cont33cso.F90  getcprj.F90                 m_gemm_nonlop.F90  nonlop_gpu.F90  opernl4b.F90     trace2.F90 \
cont33so.F90   interfaces_66_nonlocal.F90  m_hamiltonian.F90  nonlop_pl.F90   opernla_ylm.F90

include ../make.inc.post


cont13.o : cont13.F90 
cont22cso.o : cont22cso.F90 
cont22.o : cont22.F90 
cont22so.o : cont22so.F90 
cont24.o : cont24.F90 
cont33cso.o : cont33cso.F90 
cont33so.o : cont33so.F90 
cont35.o : cont35.F90 
cont3.o : cont3.F90 
ctocprj.o : ctocprj.F90 interfaces_66_nonlocal.o 
ddkten.o : ddkten.F90 
dfpt_mkffkg.o : dfpt_mkffkg.F90 
getcprj.o : getcprj.F90 interfaces_66_nonlocal.o 
interfaces_66_nonlocal.o : interfaces_66_nonlocal.F90 m_hamiltonian.o 
metcon.o : metcon.F90 
metcon_so.o : metcon_so.F90 
metric_so.o : metric_so.F90 
metstr.o : metstr.F90 
m_fock.o : m_fock.F90 
m_gemm_nonlop.o : m_gemm_nonlop.F90 interfaces_66_nonlocal.o 
m_hamiltonian.o : m_hamiltonian.F90 m_fock.o 
mkffkg.o : mkffkg.F90 
mkffnl.o : mkffnl.F90 
mkkpg.o : mkkpg.F90 
m_vkbr.o : m_vkbr.F90 interfaces_66_nonlocal.o 
nonlop.o : nonlop.F90 interfaces_66_nonlocal.o m_hamiltonian.o m_gemm_nonlop.o 
nonlop_gpu.o : nonlop_gpu.F90 
nonlop_pl.o : nonlop_pl.F90 interfaces_66_nonlocal.o 
nonlop_test.o : nonlop_test.F90 interfaces_66_nonlocal.o m_hamiltonian.o 
nonlop_ylm.o : nonlop_ylm.F90 interfaces_66_nonlocal.o 
opernl2.o : opernl2.F90 interfaces_66_nonlocal.o 
opernl3.o : opernl3.F90 interfaces_66_nonlocal.o 
opernl4a.o : opernl4a.F90 interfaces_66_nonlocal.o 
opernl4b.o : opernl4b.F90 interfaces_66_nonlocal.o 
opernla_ylm.o : opernla_ylm.F90 
opernlb_ylm.o : opernlb_ylm.F90 
opernlc_ylm.o : opernlc_ylm.F90 
opernld_ylm.o : opernld_ylm.F90 
scalewf_nonlop.o : scalewf_nonlop.F90 
strsocv.o : strsocv.F90 
trace2.o : trace2.F90 
