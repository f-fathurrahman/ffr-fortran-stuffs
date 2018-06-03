!!****m* ABINIT/interfaces_77_ddb
!! NAME
!! interfaces_77_ddb
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/77_ddb
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module interfaces_77_ddb

 implicit none

interface
 subroutine alignph(amu,displ,d2cart,mpert,natom,ntypat,phfrq,typat)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: ntypat
  real(dp),intent(in) :: amu(ntypat)
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(inout) :: displ(2,3*natom,3*natom)
  real(dp),intent(in) :: phfrq(3*natom)
  integer,intent(in) :: typat(natom)
 end subroutine alignph
end interface

interface
 subroutine complete_gamma(Cryst,nbranch,nsppol,nqptirred,nqpt_full,ep_scalprod,qirredtofull,qpttoqpt,gamma_qpt)
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: ep_scalprod
  integer,intent(in) :: nbranch
  integer,intent(in) :: nqpt_full
  integer,intent(in) :: nqptirred
  integer,intent(in) :: nsppol
  type(crystal_t),intent(in) :: Cryst
  real(dp), intent(inout) :: gamma_qpt(2,nbranch**2,nsppol,nqpt_full)
  integer,intent(in) :: qirredtofull(nqptirred)
  integer,intent(in) :: qpttoqpt(2,Cryst%nsym,nqpt_full)
 end subroutine complete_gamma
end interface

interface
 subroutine complete_gamma_tr(crystal,ep_scalprod,nbranch,nqptirred,nqpt_full,nsppol,gamma_qpt_tr,qirredtofull,qpttoqpt)
  use defs_basis
  use m_crystal
  implicit none
  integer, intent(in) :: ep_scalprod
  integer, intent(in) :: nbranch
  integer, intent(in) :: nqpt_full
  integer, intent(in) :: nqptirred
  integer, intent(in) :: nsppol
  type(crystal_t),intent(in) :: crystal
  real(dp), intent(inout) :: gamma_qpt_tr(2,9,nbranch*nbranch,nsppol,nqpt_full)
  integer,intent(in) :: qirredtofull(nqptirred)
  integer,intent(in) :: qpttoqpt(2,crystal%nsym,nqpt_full)
 end subroutine complete_gamma_tr
end interface

interface
 subroutine complete_gkk(elph_ds,gkk_flag,gprimd,indsym,natom,nsym,qpttoqpt,rprimd,symrec,symrel)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(inout) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch, &
  &         elph_ds%k_phon%my_nkpt,elph_ds%nsppol,elph_ds%nqpt_full)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: qpttoqpt(2,nsym,elph_ds%nqpt_full)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrec(3,3,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine complete_gkk
end interface

interface
 subroutine complete_gkq(Cryst,nbranch,nqptirred,nqpt_full,ikpt,nkpt,&  
  &  ep_scalprod,qirredtofull,qpttoqpt,kpttokpt,qcount,gkk_qpt_full)
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: ep_scalprod
  integer,intent(in) :: ikpt
  integer,intent(in) :: nbranch
  integer,intent(in) :: nkpt
  integer,intent(in) :: nqpt_full
  integer,intent(in) :: nqptirred
  type(crystal_t),intent(in) :: Cryst
  real(dp), intent(inout) :: gkk_qpt_full(2,nbranch**2,nqpt_full,nkpt)
  integer,intent(in) :: kpttokpt(2,Cryst%nsym,nkpt)
  integer,intent(inout) :: qcount(nkpt)
  integer,intent(in) :: qirredtofull(nqptirred)
  integer,intent(in) :: qpttoqpt(2,Cryst%nsym,nqpt_full)
 end subroutine complete_gkq
end interface

interface
 subroutine completeperts(Cryst,nbranch,nFSband,nkpt,nsppol,gkk_flag,h1_mat_el,h1_mat_el_sq,&  
  &  qpt,symq,qtimrev)
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: nFSband
  integer,intent(in) :: nbranch
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: qtimrev
  type(crystal_t),intent(in) :: Cryst
  integer,intent(inout) :: gkk_flag(nbranch,nbranch,nkpt,nsppol)
  real(dp),intent(in) :: h1_mat_el(2,nFSband**2,nbranch,nkpt,nsppol)
  real(dp),intent(out) :: h1_mat_el_sq(2,nFSband**2,nbranch**2,nkpt,nsppol)
  real(dp),intent(in) :: qpt(3)
  integer,intent(in) :: symq(4,2,Cryst%nsym)
 end subroutine completeperts
end interface

interface
 subroutine d2c_weights(elph_ds,elph_tr_ds)
  use defs_elphon
  implicit none
  type(elph_type),intent(inout) :: elph_ds
  type(elph_tr_type),intent(inout),optional :: elph_tr_ds
 end subroutine d2c_weights
end interface

interface
 subroutine d2c_wtq(elph_ds)
  use defs_elphon
  implicit none
  type(elph_type),intent(inout) :: elph_ds
 end subroutine d2c_wtq
end interface

interface
 subroutine ddb_diel(Crystal,amu,anaddb_dtset,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&  
  &  iout,lst,mpert,natom,nph2l,phfrq,comm,ncid)
  use defs_basis
  use m_anaddb_dataset
  use m_crystal
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: ncid
  integer,intent(in) :: nph2l
  type(crystal_t),intent(in) :: Crystal
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  real(dp),intent(in) :: amu(Crystal%ntypat)
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(out) :: dielt_rlx(3,3)
  real(dp),intent(inout) :: displ(2,3*natom,3*natom)
  real(dp),intent(out) :: epsinf(3,3)
  real(dp),intent(out) :: fact_oscstr(2,3,3*natom)
  real(dp),intent(in) :: lst(nph2l)
  real(dp),intent(in) :: phfrq(3*natom)
 end subroutine ddb_diel
end interface

interface
 subroutine ddb_elast(anaddb_dtset,crystal,blkval,compl,compl_clamped,compl_stress,d2asr,&  
  &  elast,elast_clamped,elast_stress,iblok,iblok_stress,&  
  &  instrain,iout,mpert,msize,natom,nblok)
  use m_anaddb_dataset
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: iblok
  integer,intent(in) :: iblok_stress
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  type(crystal_t),intent(in) :: crystal
  real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
  real(dp),intent(out) :: compl(6,6)
  real(dp),intent(out) :: compl_clamped(6,6)
  real(dp),intent(out) :: compl_stress(6,6)
  real(dp),intent(in) :: d2asr(2,3,natom,3,natom)
  real(dp),intent(out) :: elast(6,6)
  real(dp),intent(out) :: elast_clamped(6,6)
  real(dp),intent(out) :: elast_stress(6,6)
  real(dp),intent(in) :: instrain(3*natom,6)
 end subroutine ddb_elast
end interface

interface
 subroutine ddb_hybrid(acell,asr,atmfrc,dielt,dipdip,dyew,dyewq0,&  
  &  gmet,gprim,iout,natom,nrpt,rcan,rmet,&  
  &  rprim,rpt,ucvol,wghatm,xred,zeff)
  use defs_basis
  implicit none
  integer,intent(in) :: asr
  integer,intent(in) :: dipdip
  integer,intent(in) :: iout
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: acell(3)
  real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)
  real(dp),intent(inout) :: dielt(3,3)
  real(dp),intent(inout) :: dyew(2,3,natom,3,natom)
  real(dp),intent(inout) :: dyewq0(3,3,natom)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: rcan(3,natom)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(inout) :: zeff(3,3,natom)
 end subroutine ddb_hybrid
end interface

interface
 subroutine ddb_internalstr(asr,crystal,blkval,asrq0,d2asr,iblok,instrain,iout,mpert,msize,natom,nblok)
  use defs_basis
  use m_ddb
  use m_crystal
  implicit none
  integer,intent(in) :: asr
  integer,intent(in) :: iblok
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  type(asrq0_t),intent(inout) :: asrq0
  type(crystal_t),intent(in) :: crystal
  real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
  real(dp),intent(in) :: d2asr(2,3,natom,3,natom)
  real(dp),intent(out) :: instrain(3*natom,6)
 end subroutine ddb_internalstr
end interface

interface
 subroutine ddb_interpolate(ifc, crystal, inp, ddb, ddb_hdr, asrq0,&  
  &  prefix, comm)
  use m_ifc
  use m_ddb
  use m_anaddb_dataset
  use m_ddb_hdr
  use m_crystal
  implicit none
  integer,intent(in) :: comm
  type(asrq0_t),intent(inout) :: asrq0
  type(crystal_t),intent(in) :: crystal
  type(ddb_type),intent(inout) :: ddb
  type(ddb_hdr_type),intent(inout) :: ddb_hdr
  type(ifc_type),intent(in) :: ifc
  type(anaddb_dataset_type),target,intent(inout) :: inp
  character(len=*),intent(in) :: prefix
 end subroutine ddb_interpolate
end interface

interface
 subroutine ddb_piezo(anaddb_dtset,blkval,dielt_rlx,elast,iblok,instrain,iout,mpert,natom,nblok,piezo,ucvol)
  use m_anaddb_dataset
  use defs_basis
  implicit none
  integer,intent(in) :: iblok
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  integer,intent(in) :: nblok
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: blkval(2,3,mpert,3,mpert,nblok)
  real(dp),intent(in) :: dielt_rlx(3,3)
  real(dp),intent(in) :: elast(6,6)
  real(dp),intent(in) :: instrain(3*natom,6)
  real(dp),intent(out) :: piezo(6,3)
 end subroutine ddb_piezo
end interface

interface
 subroutine electrooptic(dchide,dieflag,epsinf,fact_oscstr,natom,phfrq,prtmbm,rsus,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: dieflag
  integer,intent(in) :: natom
  integer,intent(in) :: prtmbm
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: dchide(3,3,3)
  real(dp),intent(in) :: epsinf(3,3)
  real(dp),intent(in) :: fact_oscstr(2,3,3*natom)
  real(dp),intent(in) :: phfrq(3*natom)
  real(dp),intent(in) :: rsus(3*natom,3,3)
 end subroutine electrooptic
end interface

interface
 subroutine eliashberg_1d(a2f_1d,elph_ds,mustar)
  use defs_elphon
  use defs_basis
  implicit none
  type(elph_type),intent(in) :: elph_ds
  real(dp),intent(in) :: mustar
  real(dp),intent(in) :: a2f_1d(elph_ds%na2f)
 end subroutine eliashberg_1d
end interface

interface
 subroutine eli_app_m_1d (delta_1d,lambda_1d,nmatsu,z_1d)
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: z_1d(-nmatsu:nmatsu)
 end subroutine eli_app_m_1d
end interface

interface
 subroutine eli_diag_m_1d (delta_1d,lambda_1d,maxeigval,mustar,nmatsu,tc,z_1d)
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  real(dp),intent(out) :: maxeigval
  real(dp),intent(in) :: mustar
  real(dp),intent(in) :: tc
  real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: z_1d(-nmatsu:nmatsu)
 end subroutine eli_diag_m_1d
end interface

interface
 subroutine eli_lambda_1d (a2f_1d,elph_ds,lambda_1d,nmatsu,tc)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  type(elph_type),intent(in) :: elph_ds
  real(dp),intent(in) :: tc
  real(dp),intent(in) :: a2f_1d(elph_ds%na2f)
  real(dp),intent(out) :: lambda_1d(-nmatsu:nmatsu)
 end subroutine eli_lambda_1d
end interface

interface
 subroutine eli_m_iter_1d (delta_1d,lambda_1d,maxeigval,nmatsu,z_1d)
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  real(dp),intent(out) :: maxeigval
  real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
  real(dp),intent(in) :: z_1d(-nmatsu:nmatsu)
 end subroutine eli_m_iter_1d
end interface

interface
 subroutine eli_z_1d (lambda_1d,nmatsu,z_1d)
  use defs_basis
  implicit none
  integer,intent(in) :: nmatsu
  real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
  real(dp),intent(out) :: z_1d(-nmatsu:nmatsu)
 end subroutine eli_z_1d
end interface

interface
 subroutine elphon(anaddb_dtset,Cryst,Ifc,filnam,comm)
  use m_anaddb_dataset
  use defs_basis
  use m_ifc
  use m_crystal
  implicit none
  integer,intent(in) :: comm
  type(crystal_t),intent(in) :: Cryst
  type(ifc_type),intent(inout) :: Ifc
  type(anaddb_dataset_type),intent(inout) :: anaddb_dtset
  character(len=fnlen),intent(in) :: filnam(7)
 end subroutine elphon
end interface

interface
 subroutine ep_el_weights(ep_b_min, ep_b_max, eigenGS, elphsmear, enemin, enemax, nene, gprimd,&  
  &  irredtoGS, kptrlatt, max_occ, minFSband, nband, nFSband, nsppol, telphint, k_obj, tmp_wtk)
  use defs_elphon
  use defs_basis
  implicit none
  integer, intent(in) :: ep_b_max
  integer, intent(in) :: ep_b_min
  integer, intent(in) :: minFSband
  integer, intent(in) :: nFSband
  integer,intent(in) :: nband
  integer,intent(in) :: nene
  integer, intent(in) :: nsppol
  integer, intent(in) :: telphint
  real(dp), intent(in) :: elphsmear
  real(dp), intent(in) :: enemax
  real(dp), intent(in) :: enemin
  type(elph_kgrid_type), intent(in) :: k_obj
  real(dp), intent(in) :: max_occ
  integer, intent(in) :: kptrlatt(3,3)
  real(dp), intent(in) :: eigenGS(nband,k_obj%nkptirr,nsppol)
  real(dp), intent(in) :: gprimd(3,3)
  integer, intent(in) :: irredtoGS(k_obj%nkptirr)
  real(dp), intent(out) :: tmp_wtk(nFSband,k_obj%nkpt,nsppol,nene)
 end subroutine ep_el_weights
end interface

interface
 subroutine ep_fs_weights(ep_b_min, ep_b_max, eigenGS, elphsmear, fermie, gprimd,&  
  &  irredtoGS, kptrlatt, max_occ, minFSband, nband, nFSband, nsppol, telphint, k_obj)
  use defs_elphon
  use defs_basis
  implicit none
  integer, intent(in) :: ep_b_max
  integer, intent(in) :: ep_b_min
  integer, intent(in) :: minFSband
  integer, intent(in) :: nFSband
  integer,intent(in) :: nband
  integer, intent(in) :: nsppol
  integer, intent(in) :: telphint
  real(dp), intent(in) :: elphsmear
  real(dp), intent(in) :: fermie
  type(elph_kgrid_type), intent(inout) :: k_obj
  real(dp), intent(in) :: max_occ
  integer, intent(in) :: kptrlatt(3,3)
  real(dp), intent(in) :: eigenGS(nband,k_obj%nkptirr,nsppol)
  real(dp), intent(in) :: gprimd(3,3)
  integer, intent(in) :: irredtoGS(k_obj%nkptirr)
 end subroutine ep_fs_weights
end interface

interface
 subroutine ep_ph_weights(phfrq,elphsmear,omega_min,omega_max,nomega,gprimd,kptrlatt,nbranch,telphint,k_obj,tmp_wtq)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: nbranch
  integer, intent(in) :: nomega
  integer, intent(in) :: telphint
  real(dp), intent(in) :: elphsmear
  type(elph_kgrid_type), intent(inout) :: k_obj
  real(dp), intent(in) :: omega_max
  real(dp), intent(in) :: omega_min
  integer, intent(in) :: kptrlatt(3,3)
  real(dp), intent(in) :: gprimd(3,3)
  real(dp), intent(in) :: phfrq(nbranch,k_obj%nkpt)
  real(dp), intent(out) :: tmp_wtq(nbranch,k_obj%nkpt,nomega)
 end subroutine ep_ph_weights
end interface

interface
 subroutine ep_setupqpt (elph_ds,crystal,anaddb_dtset,qptrlatt,timrev)
  use defs_elphon
  use m_anaddb_dataset
  use m_crystal
  implicit none
  integer, intent(in) :: timrev
  type(anaddb_dataset_type), intent(in) :: anaddb_dtset
  type(crystal_t),intent(in) :: crystal
  type(elph_type), intent(inout) :: elph_ds
  integer, intent(out) :: qptrlatt(3,3)
 end subroutine ep_setupqpt
end interface

interface
 subroutine freeze_displ_allmodes(displ, freeze_displ, natom, outfile_radix, phfreq,&  
  &  qphon, rprimd, typat, xcart, znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  real(dp), intent(in) :: freeze_displ
  character(len=*),intent(in) :: outfile_radix
  real(dp),intent(in) :: displ(2,3*natom,3*natom)
  real(dp),intent(in) :: phfreq(3*natom)
  real(dp),intent(in) :: qphon(3)
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xcart(3,natom)
  real(dp),intent(in) :: znucl(:)
 end subroutine freeze_displ_allmodes
end interface

interface
 subroutine ftgkk (wghatm,gkk_qpt,gkk_rpt,gkqwrite,gkrwrite,gprim,ikpt_phon0,&  
  &  natom,nkpt_phon,ngkkband,nkpt_used,nqpt,nrpt,nsppol,&  
  &  qtor,rpt,qpt_full,unit_gkk_rpt,unitgkq)
  use defs_basis
  implicit none
  integer,intent(in) :: gkqwrite
  integer,intent(in) :: gkrwrite
  integer,intent(in) :: ikpt_phon0
  integer,intent(in) :: natom
  integer,intent(in) :: ngkkband
  integer,intent(in) :: nkpt_phon
  integer,intent(in) :: nkpt_used
  integer,intent(in) :: nqpt
  integer,intent(in) :: nrpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: qtor
  integer,intent(in) :: unit_gkk_rpt
  integer,intent(in) :: unitgkq
  real(dp),intent(inout) :: gkk_qpt(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_used,nsppol,nqpt)
  real(dp),intent(inout) :: gkk_rpt(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_used,nsppol,nrpt)
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qpt_full(3,nqpt)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine ftgkk
end interface

interface
 subroutine fxgkkphase(elph_ds,gkk_flag,h1_mat_el,iqptfull)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptfull
  type(elph_type),intent(in) :: elph_ds
  integer,intent(in) :: gkk_flag(elph_ds%nbranch,elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  real(dp),intent(inout) :: h1_mat_el(2,elph_ds%nFSband,elph_ds%nFSband, &
  &         elph_ds%nbranch,elph_ds%k_phon%nkpt)
 end subroutine fxgkkphase
end interface

interface
 subroutine gam_mult_displ(nbranch, displ_red, gam_bare, gam_now)
  use defs_basis
  implicit none
  integer, intent(in) :: nbranch
  real(dp), intent(in) :: displ_red(2,nbranch,nbranch)
  real(dp), intent(in) :: gam_bare(2,nbranch,nbranch)
  real(dp), intent(out) :: gam_now(2,nbranch,nbranch)
 end subroutine gam_mult_displ
end interface

interface
 subroutine get_all_gkk2(crystal,ifc,elph_ds,kptirr_phon,kpt_phon)
  use defs_elphon
  use m_ifc
  use defs_basis
  use m_crystal
  implicit none
  type(crystal_t),intent(in) :: crystal
  type(elph_type),intent(inout) :: elph_ds
  type(ifc_type),intent(in) :: ifc
  real(dp),intent(in) :: kpt_phon(3,elph_ds%k_phon%nkpt)
  real(dp),intent(in) :: kptirr_phon(3,elph_ds%k_phon%nkptirr)
 end subroutine get_all_gkk2
end interface

interface
 subroutine get_all_gkq (elph_ds,Cryst,ifc,Bst,FSfullpqtofull,nband,n1wf,onegkksize,&  
  &  qpttoqpt,ep_prt_yambo,unitgkk,ifltransport)
  use defs_elphon
  use m_ifc
  use defs_datatypes
  use m_crystal
  implicit none
  integer,intent(in) :: ep_prt_yambo
  integer,intent(in) :: ifltransport
  integer,intent(in) :: n1wf
  integer,intent(in) :: nband
  integer,intent(in) :: onegkksize
  integer,intent(in) :: unitgkk
  type(ebands_t),intent(in) :: Bst
  type(crystal_t),intent(in) :: Cryst
  type(elph_type),intent(inout) :: elph_ds
  type(ifc_type),intent(in) :: ifc
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  integer,intent(in) :: qpttoqpt(2,Cryst%nsym,elph_ds%nqpt_full)
 end subroutine get_all_gkq
end interface

interface
 subroutine get_all_gkr (elph_ds,gprim,natom,nrpt,onegkksize,rpt,qpt_full,wghatm)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  integer,intent(in) :: onegkksize
  type(elph_type),intent(inout) :: elph_ds
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qpt_full(3,elph_ds%nqpt_full)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine get_all_gkr
end interface

interface
 subroutine get_fs_bands(eigenGS,hdr,fermie,ep_b_min,ep_b_max,minFSband,maxFSband,nkptirr)
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: ep_b_max
  integer, intent(in) :: ep_b_min
  integer,intent(out) :: maxFSband
  integer,intent(out) :: minFSband
  integer,intent(out) :: nkptirr
  real(dp),intent(in) :: fermie
  type(hdr_type),intent(in) :: hdr
  real(dp),intent(in) :: eigenGS(hdr%nband(1),hdr%nkpt,hdr%nsppol)
 end subroutine get_fs_bands
end interface

interface
 subroutine get_nv_fs_en(crystal,ifc,elph_ds,eigenGS,max_occ,elph_tr_ds,omega_max)
  use defs_elphon
  use m_ifc
  use defs_basis
  use m_crystal
  implicit none
  type(crystal_t),intent(in) :: crystal
  type(elph_type),intent(inout) :: elph_ds
  type(elph_tr_type),intent(inout) :: elph_tr_ds
  type(ifc_type),intent(in) :: ifc
  real(dp), intent(in) :: max_occ
  real(dp), intent(out) :: omega_max
  real(dp), intent(in) :: eigenGS(elph_ds%nband,elph_ds%k_fine%nkptirr,elph_ds%nsppol)
 end subroutine get_nv_fs_en
end interface

interface
 subroutine get_nv_fs_temp(elph_ds,BSt,eigenGS,gprimd,max_occ,elph_tr_ds)
  use defs_elphon
  use defs_basis
  use defs_datatypes
  implicit none
  type(ebands_t),intent(inout) :: BSt
  type(elph_type),intent(inout) :: elph_ds
  type(elph_tr_type),intent(inout) :: elph_tr_ds
  real(dp), intent(in) :: max_occ
  real(dp), intent(in) :: eigenGS(elph_ds%nband,elph_ds%k_fine%nkptirr,elph_ds%nsppol)
  real(dp), intent(in) :: gprimd(3,3)
 end subroutine get_nv_fs_temp
end interface

interface
 subroutine get_tau_k(Cryst,ifc,Bst,elph_ds,elph_tr_ds,eigenGS,max_occ)
  use defs_elphon
  use m_ifc
  use defs_datatypes
  use defs_basis
  use m_crystal
  implicit none
  type(ebands_t),intent(inout) :: Bst
  type(crystal_t),intent(in) :: Cryst
  type(elph_type),intent(inout) :: elph_ds
  type(elph_tr_type), intent(inout) :: elph_tr_ds
  type(ifc_type),intent(in) :: ifc
  real(dp),intent(in) :: max_occ
  real(dp),intent(in) :: eigenGS(elph_ds%nband,elph_ds%k_phon%nkpt,elph_ds%nsppol)
 end subroutine get_tau_k
end interface

interface
 subroutine get_veloc_tr(elph_ds,elph_tr_ds)
  use defs_elphon
  implicit none
  type(elph_type),intent(in) :: elph_ds
  type(elph_tr_type),intent(inout) :: elph_tr_ds
 end subroutine get_veloc_tr
end interface

interface
 subroutine harmonic_thermo(Ifc,Crystal,amu,anaddb_dtset,iout,outfilename_radix,comm,&  
  &  thmflag)
  use m_ifc
  use m_anaddb_dataset
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: iout
  integer,intent(in),optional :: thmflag
  type(crystal_t),intent(in) :: Crystal
  type(ifc_type),intent(in) :: Ifc
  type(anaddb_dataset_type),intent(in) :: anaddb_dtset
  character(len=*),intent(in) :: outfilename_radix
  real(dp),intent(in) :: amu(Crystal%ntypat)
 end subroutine harmonic_thermo
end interface

interface
 subroutine integrate_gamma(elph_ds,FSfullpqtofull)
  use defs_elphon
  implicit none
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 end subroutine integrate_gamma
end interface

interface
 subroutine integrate_gamma_alt(elph_ds,elph_tr_ds,Cryst,gprim,kptrlatt,&  
  &  natom,nrpt,nsym,qpttoqpt,rpt,wghatm)
  use defs_elphon
  use defs_basis
  use m_crystal
  implicit none
  integer, intent(in) :: natom
  integer, intent(in) :: nrpt
  integer, intent(in) :: nsym
  type(crystal_t),intent(in) :: Cryst
  type(elph_type),intent(inout) :: elph_ds
  type(elph_tr_type), intent(inout) :: elph_tr_ds
  integer, intent(in) :: kptrlatt(3,3)
  real(dp), intent(in) :: gprim(3,3)
  integer,intent(in) :: qpttoqpt(2,nsym,elph_ds%nqpt_full)
  real(dp), intent(in) :: rpt(3,nrpt)
  real(dp), intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine integrate_gamma_alt
end interface

interface
 subroutine integrate_gamma_tr(elph_ds,FSfullpqtofull,s1,s2,&  
  &  veloc_sq1,veloc_sq2,elph_tr_ds)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: s1
  integer,intent(in) :: s2
  type(elph_type),intent(in) :: elph_ds
  type(elph_tr_type), intent(inout) :: elph_tr_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  real(dp),intent(in) :: veloc_sq1(3,elph_ds%nsppol)
  real(dp),intent(in) :: veloc_sq2(3,elph_ds%nsppol)
 end subroutine integrate_gamma_tr
end interface

interface
 subroutine integrate_gamma_tr_lova(elph_ds,FSfullpqtofull,elph_tr_ds)
  use defs_elphon
  implicit none
  type(elph_type),intent(in) :: elph_ds
  type(elph_tr_type), intent(inout) :: elph_tr_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 end subroutine integrate_gamma_tr_lova
end interface

interface
 subroutine interpolate_gkk(crystal,ifc,elph_ds,kpt_phon)
  use defs_elphon
  use m_ifc
  use defs_basis
  use m_crystal
  implicit none
  type(crystal_t),intent(in) :: crystal
  type(elph_type),intent(inout) :: elph_ds
  type(ifc_type),intent(in) :: ifc
  real(dp),intent(in) :: kpt_phon(3,elph_ds%k_phon%nkpt)
 end subroutine interpolate_gkk
end interface

interface
 subroutine k_neighbors (kpt, kptrlatt,kptrank_t, rel_kpt, kpt_phon_indices)
  use defs_basis
  use m_kptrank
  implicit none
  type(kptrank_type), intent(in) :: kptrank_t
  integer, intent(out) :: kpt_phon_indices(8)
  integer, intent(in) :: kptrlatt(3,3)
  real(dp), intent(in) :: kpt(3)
  real(dp), intent(out) :: rel_kpt(3)
 end subroutine k_neighbors
end interface

interface
 subroutine lin_interpq_gam(gamma_qpt,nbranch,nqbz,nsppol,gam_now,isppol,kptrlatt,qpt)
  use defs_basis
  implicit none
  integer, intent(in) :: isppol
  integer, intent(in) :: nbranch
  integer, intent(in) :: nqbz
  integer, intent(in) :: nsppol
  integer, intent(in) :: kptrlatt(3,3)
  real(dp), intent(out) :: gam_now(2,nbranch**2)
  real(dp),intent(in) :: gamma_qpt(2,nbranch**2,nsppol,nqbz)
  real(dp), intent(in) :: qpt(3)
 end subroutine lin_interpq_gam
end interface

interface
 subroutine mblktyp1(chkopt,ddbun,dscrpt,filnam,mddb,msym,nddb,vrsddb)
  use defs_basis
  implicit none
  integer,intent(in) :: chkopt
  integer,intent(in) :: ddbun
  integer,intent(in) :: mddb
  integer,intent(out) :: msym
  integer,intent(in) :: nddb
  integer,intent(in) :: vrsddb
  character(len=fnlen),intent(in) :: dscrpt
  character(len=fnlen),intent(in) :: filnam(mddb+1)
 end subroutine mblktyp1
end interface

interface
 subroutine mblktyp5 (chkopt,ddbun,dscrpt,filnam,mddb,msym,nddb,vrsddb)
  use defs_basis
  implicit none
  integer :: chkopt
  integer,intent(in) :: ddbun
  integer,intent(in) :: mddb
  integer,intent(out) :: msym
  integer,intent(in) :: nddb
  integer,intent(in) :: vrsddb
  character(len=fnlen),intent(in) :: dscrpt
  character(len=fnlen),intent(in) :: filnam(mddb+1)
 end subroutine mblktyp5
end interface

interface
 subroutine mk_irredpert(indsym,iqptfull,irredpert,&  
  &  natom,nbranch,nqpt,nsym,qpt,qtimrev,symq,symrel)
  implicit none
  integer,intent(in) :: iqptfull
  integer,intent(in) :: natom
  integer,intent(in) :: nbranch
  integer,intent(in) :: nqpt
  integer,intent(in) :: nsym
  integer,intent(in) :: qtimrev
  integer,intent(in) :: qpt(3)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(out) :: irredpert(7,nbranch,nbranch,nqpt)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrel(3,3,nsym)
 end subroutine mk_irredpert
end interface

interface
 subroutine mka2f(Cryst,ifc,a2f_1d,dos_phon,elph_ds,kptrlatt,mustar)
  use defs_elphon
  use m_ifc
  use defs_basis
  use m_crystal
  implicit none
  type(crystal_t),intent(in) :: Cryst
  type(elph_type),target,intent(inout) :: elph_ds
  type(ifc_type),intent(in) :: ifc
  real(dp),intent(in) :: mustar
  integer, intent(in) :: kptrlatt(3,3)
  real(dp),intent(out) :: a2f_1d(elph_ds%na2f)
  real(dp),intent(out) :: dos_phon(elph_ds%na2f)
 end subroutine mka2f
end interface

interface
 subroutine mka2fQgrid(elph_ds,fname)
  use defs_elphon
  use defs_basis
  implicit none
  type(elph_type),intent(in) :: elph_ds
  character(len=fnlen),intent(in) :: fname
 end subroutine mka2fQgrid
end interface

interface
 subroutine mka2f_tr(crystal,ifc,elph_ds,ntemper,tempermin,temperinc,pair2red,elph_tr_ds)
  use defs_elphon
  use m_ifc
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: ntemper
  type(crystal_t),intent(in) :: crystal
  type(elph_type),intent(inout) :: elph_ds
  type(elph_tr_type),intent(inout) :: elph_tr_ds
  type(ifc_type),intent(in) :: ifc
  real(dp),intent(in) :: temperinc
  real(dp),intent(in) :: tempermin
  integer,intent(in) :: pair2red(elph_ds%nenergy,elph_ds%nenergy)
 end subroutine mka2f_tr
end interface

interface
 subroutine mka2f_tr_lova(crystal,ifc,elph_ds,ntemper,tempermin,temperinc,elph_tr_ds)
  use defs_elphon
  use m_ifc
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: ntemper
  type(crystal_t),intent(in) :: crystal
  type(elph_type),intent(inout) :: elph_ds
  type(elph_tr_type),intent(inout) :: elph_tr_ds
  type(ifc_type),intent(in) :: ifc
  real(dp),intent(in) :: temperinc
  real(dp),intent(in) :: tempermin
 end subroutine mka2f_tr_lova
end interface

interface
 subroutine mkFSkgrid (elph_k, nsym, symrec, timrev)
  use defs_elphon
  implicit none
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  type(elph_kgrid_type),intent(inout) :: elph_k
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine mkFSkgrid
end interface

interface
 subroutine mkfsqgrid(kpt_phon,FStoqpt,nkpt_phon,nFSqpt,tmpFSqpt)
  use defs_basis
  implicit none
  integer,intent(out) :: nFSqpt
  integer,intent(in) :: nkpt_phon
  integer,intent(out) :: FStoqpt(nkpt_phon,nkpt_phon)
  real(dp),intent(in) :: kpt_phon(3,nkpt_phon)
  real(dp),intent(out) :: tmpFSqpt(3,nkpt_phon*nkpt_phon)
 end subroutine mkfsqgrid
end interface

interface
 subroutine mkph_linwid(Cryst,ifc,elph_ds,nqpath,qpath_vertices)
  use defs_elphon
  use m_ifc
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: nqpath
  type(crystal_t),intent(in) :: Cryst
  type(elph_type),intent(inout) :: elph_ds
  type(ifc_type),intent(in) :: ifc
  real(dp),intent(in) :: qpath_vertices(3,nqpath)
 end subroutine mkph_linwid
end interface

interface
 subroutine mkqptequiv(FSfullpqtofull,Cryst,kpt_phon,nkpt_phon,nqpt,qpttoqpt,qpt_full,mqtofull)
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: nkpt_phon
  integer,intent(in) :: nqpt
  type(crystal_t),intent(in) :: Cryst
  integer,intent(out) :: FSfullpqtofull(nkpt_phon,nqpt)
  real(dp),intent(in) :: kpt_phon(3,nkpt_phon)
  integer,intent(out),optional :: mqtofull(nqpt)
  real(dp),intent(in) :: qpt_full(3,nqpt)
  integer,intent(out) :: qpttoqpt(2,Cryst%nsym,nqpt)
 end subroutine mkqptequiv
end interface

interface
 subroutine nmsq_gam (accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&  
  &  h1_mat_el_sq,iqptirred)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptirred
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
 end subroutine nmsq_gam
end interface

interface
 subroutine nmsq_gam_sumFS(accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&  
  &  h1_mat_el_sq,iqptirred)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptirred
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
 end subroutine nmsq_gam_sumFS
end interface

interface
 subroutine nmsq_pure_gkk(accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,&  
  &  h1_mat_el_sq,iqptirred)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptirred
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
 end subroutine nmsq_pure_gkk
end interface

interface
 subroutine nmsq_pure_gkk_sumfs(accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,h1_mat_el_sq,iqptirred)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptirred
  type(elph_type),intent(in) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
 end subroutine nmsq_pure_gkk_sumfs
end interface

interface
 subroutine normsq_gkq(displ_red,eigvec,elph_ds,FSfullpqtofull,&  
  &  h1_mat_el_sq,iqptirred,phfrq_tmp,qpt_irred,qdata)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: iqptirred
  type(elph_type),intent(inout) :: elph_ds
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
  real(dp),intent(inout) :: h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband, &
  &         elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
  real(dp),intent(in) :: phfrq_tmp(elph_ds%nbranch)
  real(dp),intent(out) :: qdata(elph_ds%nbranch,elph_ds%nsppol,3)
  real(dp),intent(in) :: qpt_irred(3,elph_ds%nqptirred)
 end subroutine normsq_gkq
end interface

interface
 subroutine order_fs_kpts(kptns, nkpt, kptirr,nkptirr,FSirredtoGS)
  use defs_basis
  implicit none
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkptirr
  integer,intent(out) :: FSirredtoGS(nkptirr)
  real(dp),intent(out) :: kptirr(3,nkptirr)
  real(dp),intent(in) :: kptns(3,nkpt)
 end subroutine order_fs_kpts
end interface

interface
 subroutine outelph(elph_ds,enunit,fname)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: enunit
  type(elph_type),intent(in) :: elph_ds
  character(len=fnlen),intent(in) :: fname
 end subroutine outelph
end interface

interface
 subroutine prt_gkk_yambo(displ_cart,displ_red,kpt_phon,h1_mat_el,iqpt,&  
  &  natom,nFSband,nkpt_phon,phfrq,qptn)
  use defs_basis
  implicit none
  integer,intent(in) :: iqpt
  integer,intent(in) :: nFSband
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt_phon
  real(dp),intent(in) :: displ_cart(2,3*natom,3*natom)
  real(dp),intent(in) :: displ_red(2,3*natom,3*natom)
  real(dp),intent(in) :: h1_mat_el(2,nFSband*nFSband,3*natom,nkpt_phon,1)
  real(dp),intent(in) :: kpt_phon(3,nkpt_phon)
  real(dp),intent(in) :: phfrq(3*natom)
  real(dp),intent(in) :: qptn(3)
 end subroutine prt_gkk_yambo
end interface

interface
 subroutine ramansus(d2cart,dchide,dchidt,displ,mpert,&  
  &  natom,phfrq,qphon,qphnrm,rsus,ucvol)
  use defs_basis
  implicit none
  integer,intent(in) :: mpert
  integer,intent(in) :: natom
  real(dp),intent(in) :: qphnrm
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert)
  real(dp),intent(in) :: dchide(3,3,3)
  real(dp),intent(in) :: dchidt(natom,3,3,3)
  real(dp),intent(in) :: displ(2,3*natom,3*natom)
  real(dp),intent(in) :: phfrq(3*natom)
  real(dp),intent(inout) :: qphon(3)
  real(dp),intent(out) :: rsus(3*natom,3,3)
 end subroutine ramansus
end interface

interface
 subroutine rchkGSheader (hdr,natom,nband,unitgkk)
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(out) :: nband
  integer,intent(in) :: unitgkk
  type(hdr_type),intent(inout) :: hdr
 end subroutine rchkGSheader
end interface

interface
 subroutine read_el_veloc(nband_in,nkpt_in,kpt_in,nsppol_in,elph_tr_ds)
  use defs_elphon
  use defs_basis
  implicit none
  integer, intent(in) :: nband_in
  integer, intent(in) :: nkpt_in
  integer, intent(in) :: nsppol_in
  type(elph_tr_type), intent(inout) :: elph_tr_ds
  real(dp), intent(in) :: kpt_in(3,nkpt_in)
 end subroutine read_el_veloc
end interface

interface
 subroutine read_gkk(elph_ds,Cryst,ifc,Bst,FSfullpqtofull,gkk_flag,n1wf,nband,ep_prt_yambo,unitgkk)
  use defs_elphon
  use m_ifc
  use defs_datatypes
  use m_crystal
  implicit none
  integer,intent(in) :: ep_prt_yambo
  integer,intent(in) :: n1wf
  integer,intent(in) :: nband
  integer,intent(in) :: unitgkk
  type(ebands_t),intent(in) :: Bst
  type(crystal_t),intent(in) :: Cryst
  type(elph_type),intent(inout) :: elph_ds
  type(ifc_type),intent(in) :: ifc
  integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
  integer,intent(out) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch, &
  &         elph_ds%k_phon%my_nkpt,elph_ds%nsppol,elph_ds%nqpt_full)
 end subroutine read_gkk
end interface

interface
 subroutine relaxpol(Crystal,blkflg,blkval,etotal,fred,iatfix,iout,istrfix,&  
  &  mpert,msize,natfix,natom,nstrfix,pel,red_ptot,relaxat,relaxstr,&  
  &  strten,targetpol,usepaw)
  use defs_basis
  use m_crystal
  implicit none
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(in) :: natfix
  integer,intent(in) :: natom
  integer,intent(in) :: nstrfix
  integer,intent(in) :: relaxat
  integer,intent(in) :: relaxstr
  integer,intent(in) :: usepaw
  type(crystal_t),intent(in) :: Crystal
  real(dp),intent(in) :: etotal
  integer,intent(in) :: istrfix(6)
  integer,intent(in) :: blkflg(msize)
  real(dp),intent(inout) :: blkval(2,msize)
  real(dp),intent(in) :: fred(3,natom)
  integer,intent(in) :: iatfix(natom)
  real(dp),intent(in) :: pel(3)
  real(dp),intent(in) :: red_ptot(3)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(inout) :: targetpol(3)
 end subroutine relaxpol
end interface

interface
 subroutine symgamma(elph_ds,kphon_full2full,h1_mat_el,&  
  &  indsym,natom,nsym,symq,symrec)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  type(elph_type),intent(in) :: elph_ds
  real(dp),intent(inout) :: h1_mat_el(2,elph_ds%nFSband,elph_ds%nFSband, &
  &         elph_ds%nbranch,elph_ds%k_phon%nkpt)
  integer,intent(in) :: indsym(4,nsym,natom)
  integer,intent(in) :: kphon_full2full(2,nsym,elph_ds%k_phon%nkpt)
  integer,intent(in) :: symq(4,2,nsym)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine symgamma
end interface

interface
 subroutine test_ftgkk(elph_ds,gprim,natom,nrpt,rpt,qpt_full,wghatm)
  use defs_elphon
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nrpt
  type(elph_type),intent(inout) :: elph_ds
  real(dp),intent(in) :: gprim(3,3)
  real(dp),intent(in) :: qpt_full(3,elph_ds%nqpt_full)
  real(dp),intent(in) :: rpt(3,nrpt)
  real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 end subroutine test_ftgkk
end interface

interface
 subroutine thmeig(inp, ddb, crystal,&  
  &  elph_base_name, eig2_filnam, ddbun, iout,&  
  &  natom, mpert, msize, d2asr,&  
  &  comm)
  use m_anaddb_dataset
  use m_crystal
  use defs_basis
  use m_ddb
  implicit none
  integer,intent(in) :: comm
  integer,intent(in) :: ddbun
  integer,intent(in) :: iout
  integer,intent(in) :: mpert
  integer,intent(in) :: msize
  integer,intent(inout) :: natom
  type(crystal_t), intent(inout) :: crystal
  type(ddb_type),intent(inout) :: ddb
  character(len=*),intent(in) :: eig2_filnam
  character(len=*),intent(in) :: elph_base_name
  type(anaddb_dataset_type),intent(inout) :: inp
  real(dp),intent(inout) :: d2asr(2,3,natom,3,natom)
 end subroutine thmeig
end interface

end module interfaces_77_ddb
!!***
