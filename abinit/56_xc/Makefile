include ../make.inc

LIBNAME = lib56_xc.a

SRC = \
calc_smeared_density.F90  gammapositron.F90      m_electronpositron.F90  phase.F90           xchybrid_ncpp_cc.F90 \
dfpt_mkvxc.F90            gammapositron_fft.F90  mkcore.F90              rhohxc.F90          xcpot.F90 \
dfpt_mkvxcgga.F90         hartre.F90             m_kxc.F90               rhohxcpositron.F90 \
dfpt_mkvxc_noncoll.F90    interfaces_56_xc.F90   m_xc_vdw.F90            xcden.F90


include ../make.inc.post


calc_smeared_density.o : calc_smeared_density.F90 
dfpt_mkvxc.o : dfpt_mkvxc.F90 interfaces_56_xc.o 
dfpt_mkvxcgga.o : dfpt_mkvxcgga.F90 interfaces_56_xc.o 
dfpt_mkvxc_noncoll.o : dfpt_mkvxc_noncoll.F90 interfaces_56_xc.o 
gammapositron.o : gammapositron.F90 
gammapositron_fft.o : gammapositron_fft.F90 interfaces_56_xc.o m_electronpositron.o 
hartre.o : hartre.F90 
interfaces_56_xc.o : interfaces_56_xc.F90 m_electronpositron.o 
m_electronpositron.o : m_electronpositron.F90 
mkcore.o : mkcore.F90 
m_kxc.o : m_kxc.F90 interfaces_56_xc.o 
m_xc_vdw.o : m_xc_vdw.F90 
phase.o : phase.F90 
rhohxc.o : rhohxc.F90 interfaces_56_xc.o m_electronpositron.o m_xc_vdw.o 
rhohxcpositron.o : rhohxcpositron.F90 interfaces_56_xc.o m_electronpositron.o 
xcden.o : xcden.F90 interfaces_56_xc.o 
xchybrid_ncpp_cc.o : xchybrid_ncpp_cc.F90 interfaces_56_xc.o 
xcpot.o : xcpot.F90 interfaces_56_xc.o 
