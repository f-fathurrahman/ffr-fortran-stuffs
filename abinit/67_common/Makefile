include ../make.inc

LIBNAME = lib67_common.a

SRC = \
berryphase.F90         fresidrsp.F90             mklocl_recipspace.F90          prtefield.F90 \
berryphase_new.F90     getcgqphase.F90           mklocl_wavelets.F90            prteigrs.F90 \
calc_efg.F90           hartre1.F90               mkresi.F90                     prtene.F90 \
calc_fc.F90            initberry.F90             mkrho.F90                      prtimg.F90 \
cgq_builder.F90        init_e_field_vars.F90     mksubham.F90                   prtrhomxmn.F90 \
cgwf.F90               initmv.F90                mlwfovlp.F90                   prtxf.F90 \
clnup1.F90             initro.F90                mlwfovlp_proj.F90              prtxvf.F90 \
clnup2.F90             interfaces_67_common.F90  mlwfovlp_projpaw.F90           rhotov.F90 \
conducti_nc.F90        ionion_realspace.F90      mlwfovlp_pw.F90                scprqt.F90 \
conducti_paw_core.F90  ionion_surface.F90        mlwfovlp_radial.F90            setup1.F90 \
conducti_paw.F90       jellium.F90               mlwfovlp_seedname.F90          setup_positron.F90 \
constrf.F90            ks_ddiago.F90             mlwfovlp_setup.F90             setvtr.F90 \
dielmt.F90             linemin.F90               mlwfovlp_ylmfac.F90            spatialchempot.F90 \
dieltcel.F90           mag_constr_e.F90          mlwfovlp_ylmfar.F90            spin_current.F90 \
emispec_paw.F90        mag_constr.F90            moddiel.F90                    stress.F90 \
energy.F90             make_efg_el.F90           m_plowannier.F90               strhar.F90 \
erlxconv.F90           make_efg_ion.F90          msig.F90                       sygrad.F90 \
etheta.F90             make_grad_berry.F90       multipoles_out.F90             symrhg.F90 \
etotfor.F90            m_epjdos.F90              m_vcoul.F90                    testsusmat.F90 \
evdw_wannier.F90       m_eprenorms.F90           newkpt.F90                     uderiv.F90 \
extraprho.F90          m_ioarr.F90               nres2vres.F90                  update_e_field_vars.F90 \
extrapwf.F90           m_iowf.F90                odamix.F90                     vdw_dftd2.F90 \
fconv.F90              mkcore_inner.F90          optics_vloc.F90                vdw_dftd3.F90 \
fixsym.F90             mkcore_paw.F90            partial_dos_fractions.F90      vso_realspace_local.F90 \
forces.F90             mkcore_wvl.F90            partial_dos_fractions_paw.F90  vtorhotf.F90 \
forstr.F90             mkgrid_fft.F90            posdoppler.F90                 waveformat.F90 \
forstrnps.F90          mklocl.F90                poslifetime.F90                wvl_initro.F90 \
fresid.F90             mklocl_realspace.F90      posratecore.F90                wvl_mkrho.F90


include ../make.inc.post



berryphase.o : berryphase.F90 
berryphase_new.o : berryphase_new.F90 m_iowf.o 
calc_efg.o : calc_efg.F90 interfaces_67_common.o 
calc_fc.o : calc_fc.F90 
cgq_builder.o : cgq_builder.F90 
cgwf.o : cgwf.F90 interfaces_67_common.o 
clnup1.o : clnup1.F90 interfaces_67_common.o 
clnup2.o : clnup2.F90 
conducti_nc.o : conducti_nc.F90 interfaces_67_common.o 
conducti_paw_core.o : conducti_paw_core.F90 
conducti_paw.o : conducti_paw.F90 interfaces_67_common.o 
constrf.o : constrf.F90 interfaces_67_common.o 
dielmt.o : dielmt.F90 
dieltcel.o : dieltcel.F90 
emispec_paw.o : emispec_paw.F90 
energy.o : energy.F90 interfaces_67_common.o 
erlxconv.o : erlxconv.F90 
etheta.o : etheta.F90 
etotfor.o : etotfor.F90 interfaces_67_common.o 
evdw_wannier.o : evdw_wannier.F90 interfaces_67_common.o 
extraprho.o : extraprho.F90 interfaces_67_common.o 
extrapwf.o : extrapwf.F90 
fconv.o : fconv.F90 
fixsym.o : fixsym.F90 
forces.o : forces.F90 interfaces_67_common.o 
forstr.o : forstr.F90 interfaces_67_common.o 
forstrnps.o : forstrnps.F90 
fresid.o : fresid.F90 
fresidrsp.o : fresidrsp.F90 
getcgqphase.o : getcgqphase.F90 
hartre1.o : hartre1.F90 
initberry.o : initberry.F90 
init_e_field_vars.o : init_e_field_vars.F90 interfaces_67_common.o 
initmv.o : initmv.F90 
initro.o : initro.F90 
interfaces_67_common.o : interfaces_67_common.F90 m_epjdos.o 
ionion_realspace.o : ionion_realspace.F90 
ionion_surface.o : ionion_surface.F90 
jellium.o : jellium.F90 
ks_ddiago.o : ks_ddiago.F90 
linemin.o : linemin.F90 interfaces_67_common.o 
mag_constr_e.o : mag_constr_e.F90 
mag_constr.o : mag_constr.F90 
make_efg_el.o : make_efg_el.F90 
make_efg_ion.o : make_efg_ion.F90 
make_grad_berry.o : make_grad_berry.F90 
m_epjdos.o : m_epjdos.F90 
m_eprenorms.o : m_eprenorms.F90 
m_ioarr.o : m_ioarr.F90 
m_iowf.o : m_iowf.F90 
mkcore_inner.o : mkcore_inner.F90 
mkcore_paw.o : mkcore_paw.F90 interfaces_67_common.o 
mkcore_wvl.o : mkcore_wvl.F90 interfaces_67_common.o 
mkgrid_fft.o : mkgrid_fft.F90 
mklocl.o : mklocl.F90 interfaces_67_common.o 
mklocl_realspace.o : mklocl_realspace.F90 interfaces_67_common.o 
mklocl_recipspace.o : mklocl_recipspace.F90 
mklocl_wavelets.o : mklocl_wavelets.F90 interfaces_67_common.o 
mkresi.o : mkresi.F90 
mkrho.o : mkrho.F90 interfaces_67_common.o 
mksubham.o : mksubham.F90 
mlwfovlp.o : mlwfovlp.F90 interfaces_67_common.o 
mlwfovlp_proj.o : mlwfovlp_proj.F90 interfaces_67_common.o 
mlwfovlp_projpaw.o : mlwfovlp_projpaw.F90 interfaces_67_common.o 
mlwfovlp_pw.o : mlwfovlp_pw.F90 
mlwfovlp_radial.o : mlwfovlp_radial.F90 
mlwfovlp_seedname.o : mlwfovlp_seedname.F90 
mlwfovlp_setup.o : mlwfovlp_setup.F90 
mlwfovlp_ylmfac.o : mlwfovlp_ylmfac.F90 
mlwfovlp_ylmfar.o : mlwfovlp_ylmfar.F90 
moddiel.o : moddiel.F90 
m_plowannier.o : m_plowannier.F90 
msig.o : msig.F90 
multipoles_out.o : multipoles_out.F90 
m_vcoul.o : m_vcoul.F90 
newkpt.o : newkpt.F90 
nres2vres.o : nres2vres.F90 
odamix.o : odamix.F90 
optics_vloc.o : optics_vloc.F90 
partial_dos_fractions.o : partial_dos_fractions.F90 m_epjdos.o 
partial_dos_fractions_paw.o : partial_dos_fractions_paw.F90 m_epjdos.o 
posdoppler.o : posdoppler.F90 interfaces_67_common.o 
poslifetime.o : poslifetime.F90 
posratecore.o : posratecore.F90 
prtefield.o : prtefield.F90 
prteigrs.o : prteigrs.F90 
prtene.o : prtene.F90 
prtimg.o : prtimg.F90 interfaces_67_common.o 
prtrhomxmn.o : prtrhomxmn.F90 
prtxf.o : prtxf.F90 
prtxvf.o : prtxvf.F90 
rhotov.o : rhotov.F90 interfaces_67_common.o 
scprqt.o : scprqt.F90 interfaces_67_common.o 
setup1.o : setup1.F90 
setup_positron.o : setup_positron.F90 interfaces_67_common.o m_ioarr.o 
setvtr.o : setvtr.F90 interfaces_67_common.o 
spatialchempot.o : spatialchempot.F90 
spin_current.o : spin_current.F90 interfaces_67_common.o 
stress.o : stress.F90 interfaces_67_common.o 
strhar.o : strhar.F90 
sygrad.o : sygrad.F90 
symrhg.o : symrhg.F90 
testsusmat.o : testsusmat.F90 
uderiv.o : uderiv.F90 interfaces_67_common.o 
update_e_field_vars.o : update_e_field_vars.F90 interfaces_67_common.o 
vdw_dftd2.o : vdw_dftd2.F90 
vdw_dftd3.o : vdw_dftd3.F90 
vso_realspace_local.o : vso_realspace_local.F90 
vtorhotf.o : vtorhotf.F90 interfaces_67_common.o 
waveformat.o : waveformat.F90 
wvl_initro.o : wvl_initro.F90 
wvl_mkrho.o : wvl_mkrho.F90 
