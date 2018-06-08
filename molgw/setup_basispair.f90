!=========================================================================
subroutine setup_basispair()
  use m_eri
 implicit none
!=====
 integer :: ishell,jshell
 integer :: ibf,jbf,ijbf
!=====
  logical :: negligible_basispair

 npair = 0
 do jbf=1,nbf_eri
   do ibf=1,jbf
     if( .NOT. negligible_basispair(ibf,jbf) ) then
       npair = npair + 1
     endif
   enddo
 enddo

 call clean_allocate('index pair',index_pair_1d,(nbf_eri*(nbf_eri+1))/2)
 call clean_allocate('index basis',index_basis,2,npair)

 !
 ! Specific ordering where the first nbf pairs contain the diagonal terms ibf==jbf
 !
 npair = 0
 index_pair_1d(:) = 0
 do jbf=1,nbf_eri
   if( negligible_basispair(jbf,jbf) ) then
     call die('setup_negligible_basispair: this should not happen')
   endif
   npair = npair + 1
   index_pair_1d(jbf) = npair
   index_pair_1d(jbf) = npair
   index_basis(1,npair) = jbf
   index_basis(2,npair) = jbf
 enddo

 ijbf = nbf_eri

 do ibf=1,nbf_eri 
   do jbf=ibf+1,nbf_eri  ! Skip the diagonal terms since it is already included 
     ijbf = ijbf + 1
     if( .NOT. negligible_basispair(ibf,jbf) ) then
       npair = npair + 1
       index_pair_1d(ijbf) = npair
       index_basis(1,npair) = ibf
       index_basis(2,npair) = jbf
     endif
   enddo
 enddo


end subroutine setup_basispair
