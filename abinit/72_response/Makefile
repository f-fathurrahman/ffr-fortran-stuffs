include ../make.inc

LIBNAME = lib72_response.a

SRC = \
corrmetalwf1.F90   dfptff_die.F90        dfpt_nsteltwf.F90  eig2tot.F90                 m_eig2d.F90 \
d2frnl.F90         dfptff_ebp.F90        dfpt_nstpaw.F90    eigen_meandege.F90          m_paral_pert.F90 \
d2kindstr2.F90     dfptff_edie.F90       dfpt_nstwf.F90     elph2_fanddw.F90            m_rf2.F90 \
ddb_io_out.F90     dfptff_gbefd.F90      dfpt_prtene.F90    elt_ewald.F90               m_strain.F90 \
dfpt_accrho.F90    dfptff_gradberry.F90  dfpt_prtph.F90     eltxccore.F90               newfermie1.F90 \
dfpt_cgwf.F90      dfptff_initberry.F90  dfpt_rhofermi.F90  gaugetransfo.F90            outbsd.F90 \
dfpt_dyfro.F90     dfpt_gatherdy.F90     dfpt_rhotov.F90    getshell.F90                outddbnc.F90 \
dfpt_dyout.F90     dfpt_mkcore.F90       dfpt_scfcv.F90     hartrestr.F90               outgkk.F90 \
dfpt_dyxc1.F90     dfpt_mkrho.F90        dfpt_sydy.F90      inpgkk.F90                  qmatrix.F90 \
dfpt_eltfrhar.F90  dfpt_mkvxcstr.F90     dfpt_sygra.F90     interfaces_72_response.F90  redgr.F90 \
dfpt_eltfrkin.F90  dfpt_mkvxcstrgga.F90  dfpt_symph.F90     kptfine_av.F90              rf2_init.F90 \
dfpt_eltfrloc.F90  dfpt_newvtr.F90       dfpt_vlocal.F90    m_ddb.F90                   rotate_rho.F90 \
dfpt_eltfrxc.F90   dfptnl_doutput.F90    dfpt_vtorho.F90    m_ddb_hdr.F90               smeared_delta.F90 \
dfpt_etot.F90      dfptnl_resp.F90       dfpt_vtowfk.F90    m_ddk.F90                   vlocalstr.F90 \
dfpt_ewald.F90     dfpt_nselt.F90        dfpt_wfkfermi.F90  m_dvdb.F90                  wrtloctens.F90 \
dfptff_bec.F90     dfpt_nstdy.F90        eig2stern.F90      m_efmas.F90

include ../make.inc.post


corrmetalwf1.o : corrmetalwf1.F90 
d2frnl.o : d2frnl.F90 m_efmas.o 
d2kindstr2.o : d2kindstr2.F90 
ddb_io_out.o : ddb_io_out.F90 
dfpt_accrho.o : dfpt_accrho.F90 
dfpt_cgwf.o : dfpt_cgwf.F90 m_rf2.o 
dfpt_dyfro.o : dfpt_dyfro.F90 interfaces_72_response.o 
dfpt_dyout.o : dfpt_dyout.F90 
dfpt_dyxc1.o : dfpt_dyxc1.F90 interfaces_72_response.o 
dfpt_eltfrhar.o : dfpt_eltfrhar.F90 
dfpt_eltfrkin.o : dfpt_eltfrkin.F90 interfaces_72_response.o 
dfpt_eltfrloc.o : dfpt_eltfrloc.F90 
dfpt_eltfrxc.o : dfpt_eltfrxc.F90 interfaces_72_response.o 
dfpt_etot.o : dfpt_etot.F90 
dfpt_ewald.o : dfpt_ewald.F90 
dfptff_bec.o : dfptff_bec.F90 
dfptff_die.o : dfptff_die.F90 
dfptff_ebp.o : dfptff_ebp.F90 
dfptff_edie.o : dfptff_edie.F90 
dfptff_gbefd.o : dfptff_gbefd.F90 
dfptff_gradberry.o : dfptff_gradberry.F90 
dfptff_initberry.o : dfptff_initberry.F90 
dfpt_gatherdy.o : dfpt_gatherdy.F90 
dfpt_mkcore.o : dfpt_mkcore.F90 
dfpt_mkrho.o : dfpt_mkrho.F90 
dfpt_mkvxcstr.o : dfpt_mkvxcstr.F90 interfaces_72_response.o 
dfpt_mkvxcstrgga.o : dfpt_mkvxcstrgga.F90 
dfpt_newvtr.o : dfpt_newvtr.F90 
dfptnl_doutput.o : dfptnl_doutput.F90 m_ddb.o 
dfptnl_resp.o : dfptnl_resp.F90 
dfpt_nselt.o : dfpt_nselt.F90 interfaces_72_response.o 
dfpt_nstdy.o : dfpt_nstdy.F90 interfaces_72_response.o 
dfpt_nsteltwf.o : dfpt_nsteltwf.F90 
dfpt_nstpaw.o : dfpt_nstpaw.F90 interfaces_72_response.o 
dfpt_nstwf.o : dfpt_nstwf.F90 interfaces_72_response.o 
dfpt_prtene.o : dfpt_prtene.F90 
dfpt_prtph.o : dfpt_prtph.F90 
dfpt_rhofermi.o : dfpt_rhofermi.F90 interfaces_72_response.o 
dfpt_rhotov.o : dfpt_rhotov.F90 interfaces_72_response.o 
dfpt_scfcv.o : dfpt_scfcv.F90 interfaces_72_response.o m_rf2.o 
dfpt_sydy.o : dfpt_sydy.F90 
dfpt_sygra.o : dfpt_sygra.F90 
dfpt_symph.o : dfpt_symph.F90 
dfpt_vlocal.o : dfpt_vlocal.F90 
dfpt_vtorho.o : dfpt_vtorho.F90 interfaces_72_response.o 
dfpt_vtowfk.o : dfpt_vtowfk.F90 interfaces_72_response.o m_rf2.o 
dfpt_wfkfermi.o : dfpt_wfkfermi.F90 interfaces_72_response.o 
eig2stern.o : eig2stern.F90 interfaces_72_response.o 
eig2tot.o : eig2tot.F90 interfaces_72_response.o m_ddb_hdr.o m_ddb.o m_eig2d.o 
eigen_meandege.o : eigen_meandege.F90 
elph2_fanddw.o : elph2_fanddw.F90 
elt_ewald.o : elt_ewald.F90 
eltxccore.o : eltxccore.F90 
gaugetransfo.o : gaugetransfo.F90 
getshell.o : getshell.F90 
hartrestr.o : hartrestr.F90 
inpgkk.o : inpgkk.F90 
interfaces_72_response.o : interfaces_72_response.F90 m_rf2.o 
kptfine_av.o : kptfine_av.F90 
m_ddb.o : m_ddb.F90 m_ddb_hdr.o 
m_ddb_hdr.o : m_ddb_hdr.F90 interfaces_72_response.o 
m_ddk.o : m_ddk.F90 
m_dvdb.o : m_dvdb.F90 
m_efmas.o : m_efmas.F90 
m_eig2d.o : m_eig2d.F90 
m_paral_pert.o : m_paral_pert.F90 
m_rf2.o : m_rf2.F90 
m_strain.o : m_strain.F90 
newfermie1.o : newfermie1.F90 
outbsd.o : outbsd.F90 
outddbnc.o : outddbnc.F90 
outgkk.o : outgkk.F90 
qmatrix.o : qmatrix.F90 
redgr.o : redgr.F90 
rf2_init.o : rf2_init.F90 m_rf2.o 
rotate_rho.o : rotate_rho.F90 
smeared_delta.o : smeared_delta.F90 
vlocalstr.o : vlocalstr.F90 
wrtloctens.o : wrtloctens.F90 
