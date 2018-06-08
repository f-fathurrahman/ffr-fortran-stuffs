!=========================================================================
subroutine setup_shellpair(basis)
  use m_basis_set
  use m_eri
 implicit none

 type(basis_set),intent(in) :: basis
!=====
 integer :: ishell,jshell
 integer :: ami,amj
 integer :: ishellpair,jshellpair
!=====
  
 ishellpair = 0
 jshellpair = 0
 do jshell=1,basis%nshell
   do ishell=1,jshell 
     jshellpair = jshellpair + 1
     ! skip the identified negligible shell pairs
     if( negligible_shellpair(ishell,jshell) ) cycle
     ami = basis%shell(ishell)%am
     amj = basis%shell(jshell)%am
     ishellpair = ishellpair + 1

   enddo
 enddo
 nshellpair = ishellpair
 write(stdout,'(/,1x,a,i8,a,i8)') 'Non negligible shellpairs to be computed',nshellpair,'  over a total of',jshellpair
 write(stdout,'(1x,a,f12.4)')     'Saving (%): ', REAL(jshellpair-nshellpair,dp)/REAL(jshellpair,dp)

 allocate(index_shellpair(2,nshellpair))

 ishellpair = 0
 do jshell=1,basis%nshell
   do ishell=1,jshell 
     ! skip the identified negligible shell pairs
     if( negligible_shellpair(ishell,jshell) ) cycle
     ami = basis%shell(ishell)%am
     amj = basis%shell(jshell)%am
     ishellpair = ishellpair + 1
     ! Reverse if needed the order of the shell so to maximize the angular
     ! momentum of the first shell
     if( ami >= amj ) then
       index_shellpair(1,ishellpair) = ishell
       index_shellpair(2,ishellpair) = jshell
     else
       index_shellpair(1,ishellpair) = jshell
       index_shellpair(2,ishellpair) = ishell
     endif

   enddo
 enddo


end subroutine setup_shellpair
