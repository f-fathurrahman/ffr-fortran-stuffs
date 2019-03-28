include ../make.inc

LIBNAME = lib66_wfs.a

SRC = \
chebfi.F90           getgh1c.F90            m_bandfft_kpt.F90         prep_getghc.F90              pw_orthon.F90 \
fock2ACE.F90         getgh2c.F90            m_esymm.F90               prep_index_wavef_bandpp.F90  rayleigh_ritz.F90 \
fock_ACE_getghc.F90  getghc.F90             m_invovl.F90              prep_nonlop.F90              subdiago.F90 \
fock_getghc.F90      getghcnd.F90           multithreaded_getghc.F90  prep_sort_wavef_spin.F90     wfconv.F90 \
fxphas.F90           getgsc.F90             prep_bandfft_tabs.F90     prep_wavef_sym_do.F90 \
getdc1.F90           interfaces_66_wfs.F90  prep_fourwf.F90           prep_wavef_sym_undo.F90


include ../make.inc.post


chebfi.o : chebfi.F90 interfaces_66_wfs.o m_bandfft_kpt.o m_invovl.o 
fock2ACE.o : fock2ACE.F90 interfaces_66_wfs.o m_bandfft_kpt.o 
fock_ACE_getghc.o : fock_ACE_getghc.F90 
fock_getghc.o : fock_getghc.F90 
fxphas.o : fxphas.F90 
getdc1.o : getdc1.F90 
getgh1c.o : getgh1c.F90 
getgh2c.o : getgh2c.F90 
getghc.o : getghc.F90 interfaces_66_wfs.o m_bandfft_kpt.o 
getghcnd.o : getghcnd.F90 
getgsc.o : getgsc.F90 
interfaces_66_wfs.o : interfaces_66_wfs.F90 m_bandfft_kpt.o 
m_bandfft_kpt.o : m_bandfft_kpt.F90 
m_esymm.o : m_esymm.F90 
m_invovl.o : m_invovl.F90 interfaces_66_wfs.o m_bandfft_kpt.o 
multithreaded_getghc.o : multithreaded_getghc.F90 interfaces_66_wfs.o m_bandfft_kpt.o 
prep_bandfft_tabs.o : prep_bandfft_tabs.F90 m_bandfft_kpt.o 
prep_fourwf.o : prep_fourwf.F90 interfaces_66_wfs.o m_bandfft_kpt.o 
prep_getghc.o : prep_getghc.F90 interfaces_66_wfs.o m_bandfft_kpt.o 
prep_index_wavef_bandpp.o : prep_index_wavef_bandpp.F90 
prep_nonlop.o : prep_nonlop.F90 interfaces_66_wfs.o m_bandfft_kpt.o 
prep_sort_wavef_spin.o : prep_sort_wavef_spin.F90 
prep_wavef_sym_do.o : prep_wavef_sym_do.F90 
prep_wavef_sym_undo.o : prep_wavef_sym_undo.F90 
pw_orthon.o : pw_orthon.F90 
rayleigh_ritz.o : rayleigh_ritz.F90 interfaces_66_wfs.o 
subdiago.o : subdiago.F90 
wfconv.o : wfconv.F90 interfaces_66_wfs.o 
