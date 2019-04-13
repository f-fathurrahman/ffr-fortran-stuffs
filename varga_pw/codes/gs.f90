subroutine gram_schmidt(ik,N_states,N_k)

      use PW
      implicit none
!
!      gram schmidt orthogonalization 
!
      integer  :: ik, N_states, N_k, i, j,k
      real*8   :: norm

      do i=1,N_states
        do j=1,i-1
          wave_function_c(:,i,ik)=wave_function_c(:,i,ik)- &
&           sum(Conjg(wave_function_c(:,j,ik))*wave_function_c(:,i,ik))*wave_function_c(:,j,ik)
        enddo
        norm=0.0
        do k=1,N_k
          norm=norm+conjg(wave_function_c(k,i,ik))*wave_function_c(k,i,ik)
        enddo
        wave_function_c(:,i,ik)=wave_function_c(:,i,ik)/sqrt(norm)
      enddo

end subroutine gram_schmidt
