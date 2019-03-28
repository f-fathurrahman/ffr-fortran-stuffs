include ../make.inc

LIBNAME = lib77_ddb.a

SRC = \
alignph.F90            ep_ph_weights.F90            interpolate_gkk.F90     m_primcell_ddb_info.F90 \
complete_gamma.F90     ep_setupqpt.F90              k_neighbors.F90         m_sortph.F90 \
complete_gamma_tr.F90  freeze_displ_allmodes.F90    lin_interpq_gam.F90     nmsq_gam.F90 \
complete_gkk.F90       ftgkk.F90                    m_anaddb_dataset.F90    nmsq_gam_sumfs.F90 \
complete_gkq.F90       fxgkkphase.F90               m_atprj.F90             nmsq_pure_gkk.F90 \
completeperts.F90      gam_mult_displ.F90           m_band2eps_dataset.F90  nmsq_pure_gkk_sumfs.F90 \
d2c_weights.F90        get_all_gkk2.F90             mblktyp1.F90            normsq_gkq.F90 \
d2c_wtq.F90            get_all_gkq.F90              mblktyp5.F90            order_fs_kpts.F90 \
ddb_diel.F90           get_all_gkr.F90              m_gruneisen.F90         outelph.F90 \
ddb_elast.F90          get_fs_bands.F90             m_ifc.F90               prt_gkk_yambo.F90 \
ddb_hybrid.F90         get_nv_fs_en.F90             mka2f.F90               ramansus.F90 \
ddb_internalstr.F90    get_nv_fs_temp.F90           mka2fQgrid.F90          rchkgsheader.F90 \
ddb_interpolate.F90    get_tau_k.F90                mka2f_tr.F90            read_el_veloc.F90 \
ddb_piezo.F90          get_veloc_tr.F90             mka2f_tr_lova.F90       read_gkk.F90 \
defs_elphon.F90        harmonic_thermo.F90          mkfskgrid.F90           relaxpol.F90 \
electrooptic.F90       integrate_gamma_alt.F90      mkfsqgrid.F90           symgamma.F90 \
eliashberg_1d.F90      integrate_gamma.F90          mk_irredpert.F90        test_ftgkk.F90 \
elphon.F90             integrate_gamma_tr.F90       mkph_linwid.F90         thmeig.F90 \
ep_el_weights.F90      integrate_gamma_tr_lova.F90  mkqptequiv.F90 \
ep_fs_weights.F90      interfaces_77_ddb.F90        m_phonons.F90


include ../make.inc.post


alignph.o : alignph.F90 
complete_gamma.o : complete_gamma.F90 defs_elphon.o 
complete_gamma_tr.o : complete_gamma_tr.F90 defs_elphon.o 
complete_gkk.o : complete_gkk.F90 defs_elphon.o 
complete_gkq.o : complete_gkq.F90 defs_elphon.o 
completeperts.o : completeperts.F90 defs_elphon.o 
d2c_weights.o : d2c_weights.F90 defs_elphon.o 
d2c_wtq.o : d2c_wtq.F90 defs_elphon.o 
ddb_diel.o : ddb_diel.F90 interfaces_77_ddb.o m_anaddb_dataset.o 
ddb_elast.o : ddb_elast.F90 m_anaddb_dataset.o 
ddb_hybrid.o : ddb_hybrid.F90 
ddb_internalstr.o : ddb_internalstr.F90 
ddb_interpolate.o : ddb_interpolate.F90 m_anaddb_dataset.o m_ifc.o 
ddb_piezo.o : ddb_piezo.F90 m_anaddb_dataset.o 
defs_elphon.o : defs_elphon.F90 
electrooptic.o : electrooptic.F90 
eliashberg_1d.o : eliashberg_1d.F90 interfaces_77_ddb.o defs_elphon.o 
elphon.o : elphon.F90 interfaces_77_ddb.o m_anaddb_dataset.o m_ifc.o defs_elphon.o 
ep_el_weights.o : ep_el_weights.F90 defs_elphon.o 
ep_fs_weights.o : ep_fs_weights.F90 defs_elphon.o 
ep_ph_weights.o : ep_ph_weights.F90 defs_elphon.o 
ep_setupqpt.o : ep_setupqpt.F90 m_anaddb_dataset.o defs_elphon.o 
freeze_displ_allmodes.o : freeze_displ_allmodes.F90 
ftgkk.o : ftgkk.F90 
fxgkkphase.o : fxgkkphase.F90 defs_elphon.o 
gam_mult_displ.o : gam_mult_displ.F90 
get_all_gkk2.o : get_all_gkk2.F90 interfaces_77_ddb.o m_ifc.o defs_elphon.o 
get_all_gkq.o : get_all_gkq.F90 interfaces_77_ddb.o m_ifc.o defs_elphon.o 
get_all_gkr.o : get_all_gkr.F90 interfaces_77_ddb.o defs_elphon.o 
get_fs_bands.o : get_fs_bands.F90 
get_nv_fs_en.o : get_nv_fs_en.F90 interfaces_77_ddb.o m_ifc.o defs_elphon.o 
get_nv_fs_temp.o : get_nv_fs_temp.F90 interfaces_77_ddb.o defs_elphon.o 
get_tau_k.o : get_tau_k.F90 interfaces_77_ddb.o m_ifc.o defs_elphon.o 
get_veloc_tr.o : get_veloc_tr.F90 defs_elphon.o 
harmonic_thermo.o : harmonic_thermo.F90 m_ifc.o m_anaddb_dataset.o m_sortph.o 
integrate_gamma_alt.o : integrate_gamma_alt.F90 interfaces_77_ddb.o defs_elphon.o 
integrate_gamma.o : integrate_gamma.F90 defs_elphon.o 
integrate_gamma_tr.o : integrate_gamma_tr.F90 defs_elphon.o 
integrate_gamma_tr_lova.o : integrate_gamma_tr_lova.F90 defs_elphon.o 
interfaces_77_ddb.o : interfaces_77_ddb.F90 m_ifc.o m_anaddb_dataset.o defs_elphon.o 
interpolate_gkk.o : interpolate_gkk.F90 interfaces_77_ddb.o m_ifc.o defs_elphon.o 
k_neighbors.o : k_neighbors.F90 
lin_interpq_gam.o : lin_interpq_gam.F90 defs_elphon.o 
m_anaddb_dataset.o : m_anaddb_dataset.F90 
m_atprj.o : m_atprj.F90 
m_band2eps_dataset.o : m_band2eps_dataset.F90 
mblktyp1.o : mblktyp1.F90 
mblktyp5.o : mblktyp5.F90 
m_gruneisen.o : m_gruneisen.F90 m_anaddb_dataset.o m_ifc.o 
m_ifc.o : m_ifc.F90 
mka2f.o : mka2f.F90 interfaces_77_ddb.o m_ifc.o defs_elphon.o 
mka2fQgrid.o : mka2fQgrid.F90 defs_elphon.o 
mka2f_tr.o : mka2f_tr.F90 interfaces_77_ddb.o m_ifc.o defs_elphon.o 
mka2f_tr_lova.o : mka2f_tr_lova.F90 interfaces_77_ddb.o m_ifc.o defs_elphon.o 
mkfskgrid.o : mkfskgrid.F90 defs_elphon.o 
mkfsqgrid.o : mkfsqgrid.F90 
mk_irredpert.o : mk_irredpert.F90 
mkph_linwid.o : mkph_linwid.F90 interfaces_77_ddb.o m_ifc.o defs_elphon.o 
mkqptequiv.o : mkqptequiv.F90 defs_elphon.o 
m_phonons.o : m_phonons.F90 interfaces_77_ddb.o m_anaddb_dataset.o m_ifc.o m_sortph.o m_atprj.o 
m_primcell_ddb_info.o : m_primcell_ddb_info.F90 
m_sortph.o : m_sortph.F90 
nmsq_gam.o : nmsq_gam.F90 interfaces_77_ddb.o defs_elphon.o 
nmsq_gam_sumfs.o : nmsq_gam_sumfs.F90 interfaces_77_ddb.o defs_elphon.o 
nmsq_pure_gkk.o : nmsq_pure_gkk.F90 interfaces_77_ddb.o defs_elphon.o 
nmsq_pure_gkk_sumfs.o : nmsq_pure_gkk_sumfs.F90 interfaces_77_ddb.o defs_elphon.o 
normsq_gkq.o : normsq_gkq.F90 interfaces_77_ddb.o defs_elphon.o 
order_fs_kpts.o : order_fs_kpts.F90 defs_elphon.o 
outelph.o : outelph.F90 defs_elphon.o 
prt_gkk_yambo.o : prt_gkk_yambo.F90 
ramansus.o : ramansus.F90 
rchkgsheader.o : rchkgsheader.F90 
read_el_veloc.o : read_el_veloc.F90 defs_elphon.o 
read_gkk.o : read_gkk.F90 interfaces_77_ddb.o m_ifc.o defs_elphon.o 
relaxpol.o : relaxpol.F90 
symgamma.o : symgamma.F90 defs_elphon.o 
test_ftgkk.o : test_ftgkk.F90 interfaces_77_ddb.o defs_elphon.o 
thmeig.o : thmeig.F90 m_anaddb_dataset.o 
