!=========================================================================
subroutine distribute_auxil_basis(nbf_auxil_basis)
  USE m_scalapack
  USE m_eri
  IMPLICIT NONE 

 integer,intent(in)  :: nbf_auxil_basis
!=====
 integer :: ibf
 integer :: ibf_local
 integer :: iproc
 integer :: ilocal,iglobal
 integer :: nbf_local_iproc(0:nproc_auxil-1)
!=====

 if( parallel_buffer ) then

#if 1
   
   call set_auxil_block_size(nbf_auxil_basis/(nprow_auxil*4))

   do iproc=0,nprow_auxil-1
     nbf_local_iproc(iproc) = NUMROC(nbf_auxil_basis,MBLOCK_AUXIL,iproc,first_row,nprow_auxil)
   enddo

   nauxil_3center = nbf_local_iproc(iprow_auxil)

   allocate(ibf_auxil_g(nauxil_3center))
   do ilocal=1,nauxil_3center
     ibf_auxil_g(ilocal) = INDXL2G(ilocal,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
   enddo
   allocate(ibf_auxil_l(nbf_auxil_basis))
   allocate(iproc_ibf_auxil(nbf_auxil_basis))
   do iglobal=1,nbf_auxil_basis
     ibf_auxil_l(iglobal)     = INDXG2L(iglobal,MBLOCK_AUXIL,0,first_row,nprow_auxil)
     iproc_ibf_auxil(iglobal) = INDXG2P(iglobal,MBLOCK_AUXIL,0,first_row,nprow_auxil)
   enddo

#else

  
   allocate(iproc_ibf_auxil(nbf_auxil_basis))
  
   iproc              = nproc_auxil - 1
   nbf_local_iproc(:) = 0
   do ibf=1,nbf_auxil_basis
  
     iproc = MODULO(iproc+1,nproc_auxil)
  
     iproc_ibf_auxil(ibf) = iproc
  
     nbf_local_iproc(iproc) = nbf_local_iproc(iproc) + 1
  
   enddo
  
   nauxil_3center = nbf_local_iproc(rank_auxil)
  
   allocate(ibf_auxil_g(nauxil_3center))
   allocate(ibf_auxil_l(nbf_auxil_basis))
   ibf_auxil_l(:) = 0
   ibf_local = 0
   do ibf=1,nbf_auxil_basis
     if( rank_auxil == iproc_ibf_auxil(ibf) ) then
       ibf_local = ibf_local + 1
       ibf_auxil_g(ibf_local) = ibf
       ibf_auxil_l(ibf)       = ibf_local
     endif
   enddo
#endif
  
 else

   ! Use SCALAPACK routines to distribute the auxiliary basis
   ! Assume a processor grid: nproc_auxil x 1

#if 1
   
   call set_auxil_block_size(nbf_auxil_basis/nprow_auxil/2)

   nauxil_3center = NUMROC(nbf_auxil_basis,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
   allocate(ibf_auxil_g(nauxil_3center))
   do ilocal=1,nauxil_3center
     ibf_auxil_g(ilocal) = INDXL2G(ilocal,MBLOCK_AUXIL,iprow_auxil,first_row,nprow_auxil)
   enddo
   allocate(ibf_auxil_l(nbf_auxil_basis))
   allocate(iproc_ibf_auxil(nbf_auxil_basis))
   do iglobal=1,nbf_auxil_basis
     ibf_auxil_l(iglobal)     = INDXG2L(iglobal,MBLOCK_AUXIL,0,first_row,nprow_auxil)
     iproc_ibf_auxil(iglobal) = INDXG2P(iglobal,MBLOCK_AUXIL,0,first_row,nprow_auxil)
   enddo

#else

   allocate(iproc_ibf_auxil(nbf_auxil_basis))
  
   iproc              = nproc_local-1
   nbf_local_iproc(:) = 0
   do ibf=1,nbf_auxil_basis
  
     iproc = MODULO(iproc+1,nproc_local)
  
     iproc_ibf_auxil(ibf) = iproc
  
     nbf_local_iproc(iproc) = nbf_local_iproc(iproc) + 1
  
   enddo
  
   nauxil_3center = nbf_local_iproc(rank_local)
  
   allocate(ibf_auxil_g(nauxil_3center))
   allocate(ibf_auxil_l(nbf_auxil_basis))
   ibf_auxil_l(:) = 0
   ibf_local = 0
   do ibf=1,nbf_auxil_basis
     if( rank_local == iproc_ibf_auxil(ibf) ) then
       ibf_local = ibf_local + 1
       ibf_auxil_g(ibf_local) = ibf
       ibf_auxil_l(ibf)       = ibf_local
     endif
   enddo

#endif
  
 endif

 write(stdout,'(/,a)') ' Distribute auxiliary basis functions among processors'
 write(stdout,'(1x,a,i4,a,i6,a)') 'Max auxiliary basis functions ',MAXVAL(nbf_local_iproc(:)),' for processor ',MAXLOC(nbf_local_iproc,DIM=1)
 write(stdout,'(1x,a,i4,a,i6,a)') 'Min auxiliary basis functions ',MINVAL(nbf_local_iproc(:)),' for processor ',MINLOC(nbf_local_iproc,DIM=1)


end subroutine distribute_auxil_basis
