include ../make.inc

LIBNAME = lib43_wvl_wrappers.a

SRC = \
defs_wvltypes.F90 \
defs_PSolver.F90   paw2wvl_ij.F90               wvl_descr_free.F90       wvl_setngfft.F90 \
wvl_cprjreorder.F90   wvl_descr_psp_set.F90    wvl_wfs_free.F90 \
interfaces_43_wvl_wrappers.F90  wvl_denspot_free.F90         wvl_paw_free.F90         wvl_wfs_lr_copy.F90 \
m_abi2big.F90                   wvl_denspot_set.F90          wvl_projectors_free.F90  wvl_wfs_set.F90 \
nullify_wvl_data.F90            wvl_descr_atoms_set.F90      wvl_projectors_set.F90 \
paw2wvl.F90                     wvl_descr_atoms_set_sym.F90  wvl_setboxgeometry.F90


include ../make.inc.post


defs_PSolver.o : defs_PSolver.F90 
defs_wvltypes.o : defs_wvltypes.F90 
interfaces_43_wvl_wrappers.o : interfaces_43_wvl_wrappers.F90 defs_wvltypes.o 
m_abi2big.o : m_abi2big.F90 defs_wvltypes.o 
nullify_wvl_data.o : nullify_wvl_data.F90 defs_wvltypes.o 
paw2wvl.o : paw2wvl.F90 defs_wvltypes.o 
paw2wvl_ij.o : paw2wvl_ij.F90 defs_wvltypes.o 
wvl_cprjreorder.o : wvl_cprjreorder.F90 defs_wvltypes.o 
wvl_denspot_free.o : wvl_denspot_free.F90 defs_wvltypes.o 
wvl_denspot_set.o : wvl_denspot_set.F90 defs_wvltypes.o 
wvl_descr_atoms_set.o : wvl_descr_atoms_set.F90 defs_wvltypes.o 
wvl_descr_atoms_set_sym.o : wvl_descr_atoms_set_sym.F90 defs_wvltypes.o 
wvl_descr_free.o : wvl_descr_free.F90 defs_wvltypes.o 
wvl_descr_psp_set.o : wvl_descr_psp_set.F90 defs_wvltypes.o 
wvl_paw_free.o : wvl_paw_free.F90 defs_wvltypes.o 
wvl_projectors_free.o : wvl_projectors_free.F90 defs_wvltypes.o 
wvl_projectors_set.o : wvl_projectors_set.F90 defs_wvltypes.o 
wvl_setboxgeometry.o : wvl_setboxgeometry.F90 defs_wvltypes.o 
wvl_setngfft.o : wvl_setngfft.F90 
wvl_wfs_free.o : wvl_wfs_free.F90 defs_wvltypes.o 
wvl_wfs_lr_copy.o : wvl_wfs_lr_copy.F90 defs_wvltypes.o 
wvl_wfs_set.o : wvl_wfs_set.F90 defs_wvltypes.o 
