!=========================================================================
subroutine distribute_auxil_basis_lr(nbf_auxil_basis)
 use m_scalapack
 USE m_eri
 implicit NONE

 integer,intent(in)  :: nbf_auxil_basis
!=====
 integer :: ibf
 integer :: ibf_local
 integer :: iproc
 integer :: ilocal,iglobal
 integer :: nbf_local_iproc_lr(0:nproc_auxil-1)
!=====

   
#ifdef HAVE_SCALAPACK

 do iproc=0,nprow_auxil-1
   nbf_local_iproc_lr(iproc) = NUMROC(nbf_auxil_basis,MBLOCK_AUXIL,iproc,first_row,nprow_auxil)
 enddo

 nauxil_3center_lr = nbf_local_iproc_lr(iprow_auxil)

 allocate(ibf_auxil_g_lr(nauxil_3center_lr))
 do ilocal=1,nauxil_3center_lr
   ibf_auxil_g_lr(ilocal) = INDXL2G(ilocal,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
 enddo
 allocate(ibf_auxil_l_lr(nbf_auxil_basis))
 allocate(iproc_ibf_auxil_lr(nbf_auxil_basis))
 do iglobal=1,nbf_auxil_basis
   ibf_auxil_l_lr(iglobal)     = INDXG2L(iglobal,MBLOCK_AUXIL,0,first_row,nprow_auxil)
   iproc_ibf_auxil_lr(iglobal) = INDXG2P(iglobal,MBLOCK_AUXIL,0,first_row,nprow_auxil)
 enddo

#else

 allocate(iproc_ibf_auxil_lr(nbf_auxil_basis))

 iproc = nproc_auxil - 1
 nbf_local_iproc_lr(:) = 0
 do ibf=1,nbf_auxil_basis

   iproc = MODULO(iproc+1,nproc_auxil)

   iproc_ibf_auxil_lr(ibf) = iproc

   nbf_local_iproc_lr(iproc) = nbf_local_iproc_lr(iproc) + 1

 enddo

 nauxil_3center_lr = nbf_local_iproc_lr(rank_auxil)

 allocate(ibf_auxil_g_lr(nauxil_3center_lr))
 allocate(ibf_auxil_l_lr(nbf_auxil_basis))
 ibf_auxil_l_lr(:) = 0
 ibf_local = 0
 do ibf=1,nbf_auxil_basis
   if( rank_auxil == iproc_ibf_auxil_lr(ibf) ) then
     ibf_local = ibf_local + 1
     ibf_auxil_g_lr(ibf_local) = ibf
     ibf_auxil_l_lr(ibf)       = ibf_local
   endif
 enddo

#endif

 write(stdout,'(/,a)') ' Distribute LR auxiliary basis functions among processors'
 write(stdout,'(1x,a,i4,a,i6,a)')   'Max auxiliary basis functions ',MAXVAL(nbf_local_iproc_lr(:)),' for processor ',MAXLOC(nbf_local_iproc_lr,DIM=1)
 write(stdout,'(1x,a,i4,a,i6,a)')   'Min auxiliary basis functions ',MINVAL(nbf_local_iproc_lr(:)),' for processor ',MINLOC(nbf_local_iproc_lr,DIM=1)


end subroutine distribute_auxil_basis_lr
