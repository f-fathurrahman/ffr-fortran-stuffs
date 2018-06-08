!=========================================================================
subroutine evaluate_s2_operator(nbf,nstate,occupation,c_matrix,s_matrix)
  use m_definitions
  use m_inputparam, only : nspin
 implicit none
 integer,intent(in)      :: nbf,nstate
 real(dp),intent(in)     :: occupation(nstate,nspin)
 real(dp),intent(in)     :: c_matrix(nbf,nstate,nspin)
 real(dp),intent(in)     :: s_matrix(nbf,nbf)
!=====
 integer                 :: ispin,istate,jstate
 real(dp)                :: s2,s2_exact
 real(dp)                :: n1,n2,nmax,nmin
!=====

 if(nspin /= 2) return

 n1 = SUM( occupation(:,1) )
 n2 = SUM( occupation(:,2) )
 nmax = MAX(n1,n2)
 nmin = MIN(n1,n2)

 s2_exact = (nmax-nmin)/2.0_dp * ( (nmax-nmin)/2.0_dp + 1.0_dp )
 s2       = s2_exact + nmin
 do istate=1,nstate
   if( occupation(istate,1) < completely_empty ) cycle
   do jstate=1,nstate
     if( occupation(jstate,2) < completely_empty ) cycle

     s2 = s2 - ABS( DOT_PRODUCT( c_matrix(:,istate,1) , MATMUL( s_matrix(:,:) , c_matrix(:,jstate,2) ) )  &
                      * occupation(istate,1) * occupation(jstate,2) )**2

   enddo
 enddo


 write(stdout,'(/,a,f8.4)') ' Total Spin S**2: ',s2
 write(stdout,'(a,f8.4)')   ' Instead of:      ',s2_exact


end subroutine evaluate_s2_operator
