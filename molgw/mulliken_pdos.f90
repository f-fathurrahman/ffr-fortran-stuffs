!=========================================================================
subroutine mulliken_pdos(nstate,basis,s_matrix,c_matrix,occupation,energy)
 use m_definitions
 use m_mpi
 use m_inputparam, only: nspin
 use m_atoms
 use m_basis_set
 implicit none
 integer,intent(in)         :: nstate
 type(basis_set),intent(in) :: basis
 real(dp),intent(in)        :: s_matrix(basis%nbf,basis%nbf)
 real(dp),intent(in)        :: c_matrix(basis%nbf,nstate,nspin)
 real(dp),intent(in)        :: occupation(nstate,nspin),energy(nstate,nspin)
!=====
 integer                    :: ibf,li,ibf1,ibf2,ishell
 integer                    :: natom1,natom2,istate,ispin
 logical                    :: file_exists
 integer                    :: pdosfile
 real(dp)                   :: proj_state_i(0:basis%ammax),proj_charge
 real(dp)                   :: cs_vector_i(basis%nbf)
 integer                    :: iatom_ibf(basis%nbf)
 integer                    :: li_ibf(basis%nbf)
!=====

 write(stdout,*)
 write(stdout,*) 'Projecting wavefunctions on selected atoms'
 inquire(file='manual_pdos',exist=file_exists)
 if(file_exists) then
   write(stdout,*) 'Opening file:','manual_pdos'
   open(newunit=pdosfile,file='manual_pdos',status='old')
   read(pdosfile,*) natom1,natom2
   close(pdosfile)
 else
   natom1=1
   natom2=1
 endif
 write(stdout,'(1x,a,i5,2x,i5)') 'Range of atoms considered: ',natom1,natom2

 do ishell=1,basis%nshell
   ibf1    = basis%shell(ishell)%istart
   ibf2    = basis%shell(ishell)%iend

   iatom_ibf(ibf1:ibf2) = basis%shell(ishell)%iatom
   li_ibf(ibf1:ibf2)    = basis%shell(ishell)%am
 enddo
 

 write(stdout,'(1x,a)') '==========================================='
 write(stdout,'(1x,a)') 'spin state  energy(eV)  Mulliken proj. total        proj s         proj p      proj d ... '
 proj_charge = 0.0_dp
 do ispin=1,nspin
   do istate=1,nstate
     proj_state_i(:) = 0.0_dp

     cs_vector_i(:) = MATMUL( c_matrix(:,istate,ispin) , s_matrix(:,:) )

     do ibf=1,basis%nbf
       if( iatom_ibf(ibf) >= natom1 .AND. iatom_ibf(ibf) <= natom2 ) then
         li = li_ibf(ibf)
         proj_state_i(li) = proj_state_i(li) + c_matrix(ibf,istate,ispin) * cs_vector_i(ibf)
       endif
     enddo
     proj_charge = proj_charge + occupation(istate,ispin) * SUM(proj_state_i(:))

     write(stdout,'(i3,1x,i5,1x,20(f16.6,4x))') ispin,istate,energy(istate,ispin) * Ha_eV,&
          SUM(proj_state_i(:)),proj_state_i(:)
   enddo
 enddo
 write(stdout,'(1x,a)') '==========================================='
 write(stdout,'(1x,a,f12.6)') 'Total Mulliken charge: ',proj_charge


end subroutine mulliken_pdos
