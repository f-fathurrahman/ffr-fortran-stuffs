include ../make.inc

LIBNAME = lib70_gw.a

SRC = \
calc_rpa_functional.F90  dotproductqrc.F90                  gwls_Projected_AT.F90          m_qparticles.F90 \
calc_sigc_cd.F90         gwcompleteness.F90                 gwls_Projected_BT.F90          m_screen.F90 \
calc_sigc_me.F90         gwls_ComputeCorrelationEnergy.F90  gwls_QR_factorization.F90      m_screening.F90 \
calc_sigc_pole_cd.F90    gwls_ComputePoles.F90              gwls_TimingLog.F90             m_sigma.F90 \
calc_sig_ppm_comp.F90    gwls_DielectricArray.F90           gwls_valenceWavefunctions.F90  my_calc_wfwfg.F90 \
calc_sigx_me.F90         gwls_GenerateEpsilon.F90           gwls_wf.F90                    paw_qpscgw.F90 \
calc_ucrpa.F90           gwls_GWanalyticPart.F90            gw_tools.F90                   prep_calc_ucrpa.F90 \
calc_vhxc_me.F90         gwls_GWlanczos.F90                 interfaces_70_gw.F90           random_stopping_power.F90 \
calc_wf_qp.F90           gwls_hamiltonian.F90               m_chi0.F90                     read_plowannier.F90 \
cchi0.F90                gwls_LanczosBasis.F90              m_dyson_solver.F90             setup_screening.F90 \
cchi0q0.F90              gwls_LanczosResolvents.F90         m_fft_prof.F90                 setup_sigma.F90 \
cchi0q0_intraband.F90    gwls_lineqsolver.F90               mlwfovlp_qp.F90 \
cohsex_me.F90            gwls_model_polarisability.F90      m_model_screening.F90 \
debug_tools.F90          gwls_polarisability.F90            m_ppmodel.F90


include ../make.inc.post


calc_rpa_functional.o : calc_rpa_functional.F90 
calc_sigc_cd.o : calc_sigc_cd.F90 
calc_sigc_me.o : calc_sigc_me.F90 interfaces_70_gw.o m_sigma.o m_ppmodel.o m_screening.o 
calc_sigc_pole_cd.o : calc_sigc_pole_cd.F90 m_model_screening.o 
calc_sig_ppm_comp.o : calc_sig_ppm_comp.F90 
calc_sigx_me.o : calc_sigx_me.F90 interfaces_70_gw.o m_sigma.o 
calc_ucrpa.o : calc_ucrpa.F90 
calc_vhxc_me.o : calc_vhxc_me.F90 
calc_wf_qp.o : calc_wf_qp.F90 
cchi0.o : cchi0.F90 interfaces_70_gw.o m_chi0.o 
cchi0q0.o : cchi0q0.F90 interfaces_70_gw.o m_chi0.o 
cchi0q0_intraband.o : cchi0q0_intraband.F90 m_chi0.o 
cohsex_me.o : cohsex_me.F90 interfaces_70_gw.o m_sigma.o m_screening.o 
debug_tools.o : debug_tools.F90 
dotproductqrc.o : dotproductqrc.F90 
gwcompleteness.o : gwcompleteness.F90 
gwls_ComputeCorrelationEnergy.o : gwls_ComputeCorrelationEnergy.F90 gwls_LanczosBasis.o gwls_TimingLog.o gwls_GWanalyticPart.o gwls_GenerateEpsilon.o gwls_GWlanczos.o gwls_Projected_BT.o gwls_Projected_AT.o gwls_ComputePoles.o gwls_DielectricArray.o gwls_model_polarisability.o gwls_polarisability.o gwls_lineqsolver.o gwls_hamiltonian.o gwls_wf.o 
gwls_ComputePoles.o : gwls_ComputePoles.F90 gwls_LanczosBasis.o gwls_TimingLog.o gwls_GWanalyticPart.o gwls_GenerateEpsilon.o gwls_GWlanczos.o gwls_polarisability.o gwls_lineqsolver.o gwls_hamiltonian.o gwls_wf.o 
gwls_DielectricArray.o : gwls_DielectricArray.F90 gwls_LanczosBasis.o gwls_QR_factorization.o gwls_TimingLog.o gwls_GenerateEpsilon.o gwls_model_polarisability.o gwls_polarisability.o gwls_lineqsolver.o gwls_hamiltonian.o gwls_valenceWavefunctions.o gwls_wf.o 
gwls_GenerateEpsilon.o : gwls_GenerateEpsilon.F90 gwls_GWlanczos.o gwls_model_polarisability.o gwls_polarisability.o gwls_TimingLog.o gwls_lineqsolver.o gwls_hamiltonian.o gwls_wf.o 
gwls_GWanalyticPart.o : gwls_GWanalyticPart.F90 gwls_TimingLog.o gwls_polarisability.o gwls_model_polarisability.o gwls_LanczosBasis.o gwls_GenerateEpsilon.o gwls_GWlanczos.o gwls_lineqsolver.o gwls_hamiltonian.o gwls_wf.o 
gwls_GWlanczos.o : gwls_GWlanczos.F90 gwls_QR_factorization.o gwls_polarisability.o gwls_lineqsolver.o gwls_hamiltonian.o gwls_wf.o gwls_TimingLog.o 
gwls_hamiltonian.o : gwls_hamiltonian.F90 gwls_wf.o 
gwls_LanczosBasis.o : gwls_LanczosBasis.F90 gwls_GenerateEpsilon.o gwls_hamiltonian.o gwls_wf.o 
gwls_LanczosResolvents.o : gwls_LanczosResolvents.F90 gwls_GWlanczos.o gwls_lineqsolver.o gwls_hamiltonian.o gwls_TimingLog.o gwls_wf.o 
gwls_lineqsolver.o : gwls_lineqsolver.F90 gwls_hamiltonian.o gwls_wf.o 
gwls_model_polarisability.o : gwls_model_polarisability.F90 gwls_lineqsolver.o gwls_hamiltonian.o gwls_valenceWavefunctions.o gwls_wf.o 
gwls_polarisability.o : gwls_polarisability.F90 gwls_lineqsolver.o gwls_hamiltonian.o gwls_valenceWavefunctions.o gwls_wf.o 
gwls_Projected_AT.o : gwls_Projected_AT.F90 gwls_GWanalyticPart.o gwls_LanczosResolvents.o gwls_LanczosBasis.o gwls_GWlanczos.o gwls_lineqsolver.o gwls_hamiltonian.o gwls_TimingLog.o gwls_wf.o 
gwls_Projected_BT.o : gwls_Projected_BT.F90 gwls_LanczosResolvents.o gwls_DielectricArray.o gwls_LanczosBasis.o gwls_GWlanczos.o gwls_lineqsolver.o gwls_hamiltonian.o gwls_TimingLog.o gwls_wf.o 
gwls_QR_factorization.o : gwls_QR_factorization.F90 gwls_hamiltonian.o gwls_wf.o gwls_TimingLog.o 
gwls_TimingLog.o : gwls_TimingLog.F90 gwls_hamiltonian.o 
gwls_valenceWavefunctions.o : gwls_valenceWavefunctions.F90 gwls_hamiltonian.o 
gwls_wf.o : gwls_wf.F90 
gw_tools.o : gw_tools.F90 
interfaces_70_gw.o : interfaces_70_gw.F90 m_screening.o m_ppmodel.o m_sigma.o 
m_chi0.o : m_chi0.F90 interfaces_70_gw.o 
m_dyson_solver.o : m_dyson_solver.F90 m_sigma.o 
m_fft_prof.o : m_fft_prof.F90 
mlwfovlp_qp.o : mlwfovlp_qp.F90 interfaces_70_gw.o m_qparticles.o 
m_model_screening.o : m_model_screening.F90 
m_ppmodel.o : m_ppmodel.F90 
m_qparticles.o : m_qparticles.F90 m_sigma.o 
m_screen.o : m_screen.F90 interfaces_70_gw.o m_ppmodel.o m_screening.o 
m_screening.o : m_screening.F90 
m_sigma.o : m_sigma.F90 m_screening.o 
my_calc_wfwfg.o : my_calc_wfwfg.F90 
paw_qpscgw.o : paw_qpscgw.F90 
prep_calc_ucrpa.o : prep_calc_ucrpa.F90 interfaces_70_gw.o 
random_stopping_power.o : random_stopping_power.F90 
read_plowannier.o : read_plowannier.F90 
setup_screening.o : setup_screening.F90 
setup_sigma.o : setup_sigma.F90 interfaces_70_gw.o m_screening.o
 
