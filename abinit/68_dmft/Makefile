include ../make.inc

LIBNAME = lib68_dmft.a

SRC = \
compute_levels.F90  hybridization_asymptotic_coefficient.F90  m_entropyDMFT.F90  newton.F90 \
datafordmft.F90     impurity_solve.F90                        m_green.F90        psichi_renormalization.F90 \
dmft_solve.F90      interfaces_68_dmft.F90                    m_hu.F90           qmc_prep_ctqmc.F90 \
dyson.F90           ldau_self.F90                             m_matlu.F90        spectral_function.F90 \
fermi_green.F90     local_ks_green.F90                        m_oper.F90         testcode_ctqmc.F90 \
hubbard_one.F90     m_energy.F90                              m_self.F90


include ../make.inc.post


compute_levels.o : compute_levels.F90 m_matlu.o m_oper.o 
datafordmft.o : datafordmft.F90 interfaces_68_dmft.o m_oper.o m_matlu.o 
dmft_solve.o : dmft_solve.F90 interfaces_68_dmft.o m_matlu.o m_energy.o m_hu.o m_self.o m_oper.o m_green.o 
dyson.o : dyson.F90 m_self.o m_matlu.o m_oper.o m_green.o 
fermi_green.o : fermi_green.F90 interfaces_68_dmft.o m_self.o m_green.o 
hubbard_one.o : hubbard_one.F90 interfaces_68_dmft.o m_hu.o m_matlu.o m_oper.o m_green.o 
hybridization_asymptotic_coefficient.o : hybridization_asymptotic_coefficient.F90 m_matlu.o m_oper.o 
impurity_solve.o : impurity_solve.F90 interfaces_68_dmft.o m_energy.o m_self.o m_matlu.o m_hu.o m_oper.o m_green.o 
interfaces_68_dmft.o : interfaces_68_dmft.F90 m_matlu.o m_hu.o m_green.o m_self.o m_oper.o 
ldau_self.o : ldau_self.F90 m_energy.o m_oper.o m_self.o m_green.o 
local_ks_green.o : local_ks_green.F90 m_green.o 
m_energy.o : m_energy.F90 m_oper.o m_matlu.o m_self.o m_green.o 
m_entropyDMFT.o : m_entropyDMFT.F90 
m_green.o : m_green.F90 m_self.o m_matlu.o m_oper.o 
m_hu.o : m_hu.F90 
m_matlu.o : m_matlu.F90 
m_oper.o : m_oper.F90 m_matlu.o 
m_self.o : m_self.F90 m_hu.o m_matlu.o m_oper.o 
newton.o : newton.F90 m_self.o m_green.o 
psichi_renormalization.o : psichi_renormalization.F90 m_matlu.o m_oper.o m_green.o 
qmc_prep_ctqmc.o : qmc_prep_ctqmc.F90 interfaces_68_dmft.o m_hu.o m_matlu.o m_self.o m_oper.o m_green.o 
spectral_function.o : spectral_function.F90 interfaces_68_dmft.o m_energy.o m_self.o m_hu.o m_matlu.o m_green.o 
testcode_ctqmc.o : testcode_ctqmc.F90 
