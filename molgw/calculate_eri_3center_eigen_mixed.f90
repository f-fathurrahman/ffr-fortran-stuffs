!=================================================================
subroutine calculate_eri_3center_eigen_mixed(nbf,nstate,c_matrix)
  use m_warning
  USE m_definitions
  USE m_inputparam,only: nspin
  USE m_eri
  USE m_timing
  USE m_eri_ao_mo
 implicit none

 integer,intent(in)   :: nbf,nstate
 real(dp),intent(in)  :: c_matrix(nbf,nstate,nspin)
!=====
 integer              :: kbf,lbf
 integer              :: lstate
 integer              :: klspin
 real(dp),allocatable :: eri_3center_tmp(:,:,:)
 real(dp),allocatable :: c_matrix_exx(:,:,:)
 logical              :: file_exists
!=====
  INTEGER :: index_pair
  LOGICAL :: negligible_basispair


 call start_clock(timing_eri_3center_eigen)

 inquire(file='fort.1000',exist=file_exists)
 if( .NOT. file_exists ) call die('fort.1000 not found')

 allocate(c_matrix_exx(nbf,nstate,nspin))
 open(1000,form='unformatted')
 do klspin=1,nspin
   do lstate=1,nstate
     read(1000) c_matrix_exx(:,lstate,klspin)
   enddo
 enddo
 close(1000,status='delete')


 write(stdout,'(/,a)') ' Calculate 3-center integrals on MIXED eigenstates'


 !TODO merge the 2 last indexes to save a factor 2! (i<->j symmetry)
 call clean_allocate('3-center mixed MO integrals',eri_3center_eigen_mixed,nauxil_3center,nstate,nstate,nspin)
 eri_3center_eigen_mixed(:,:,:,:) = 0.0_dp

 allocate(eri_3center_tmp(nauxil_3center,nbf,nstate)) 

 !TODO fix all this mess here to make it more similar to the previous subroutine
 do klspin=1,nspin
   ! Transformation of the first index
   eri_3center_tmp(:,:,:) = 0.0_dp
   do kbf=1,nbf
     do lbf=1,nbf
       if( negligible_basispair(kbf,lbf) ) cycle

         do lstate=1,nstate
           eri_3center_tmp(:,kbf,lstate) = eri_3center_tmp(:,kbf,lstate) &
                                      + c_matrix_exx(lbf,lstate,klspin) * eri_3center(:,index_pair(kbf,lbf))
         enddo

     enddo
   enddo
   ! Transformation of the second index
   do lstate=1,nstate
     eri_3center_eigen_mixed(:,:,lstate,klspin) = MATMUL( eri_3center_tmp(:,:,lstate) , c_matrix(:,:,klspin) )
   enddo

 enddo ! klspin
 deallocate(eri_3center_tmp)

 call stop_clock(timing_eri_3center_eigen)

end subroutine calculate_eri_3center_eigen_mixed
