! from rosetta code

subroutine Merge(A,NA,B,NB,C,NC)
 
   integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
   integer, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
   integer, intent(in)     :: B(NB)
   integer, intent(in out) :: C(NC)
 
   integer :: I,J,K
 
   I = 1; J = 1; K = 1;
   do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
         C(K) = A(I)
         I = I+1
      else
         C(K) = B(J)
         J = J+1
      endif
      K = K + 1
   enddo
   do while (I <= NA)
      C(K) = A(I)
      I = I + 1
      K = K + 1
   enddo
   return
 
end subroutine merge
 
recursive subroutine MergeSort(A,N,T)
 
   integer, intent(in) :: N
   integer, dimension(N), intent(in out) :: A
   integer, dimension((N+1)/2), intent (out) :: T
 
   integer :: NA,NB,V
 
   if (N < 2) return
   if (N == 2) then
      if (A(1) > A(2)) then
         V = A(1)
         A(1) = A(2)
         A(2) = V
      endif
      return
   endif      
   NA=(N+1)/2
   NB=N-NA
 
   call MergeSort(A,NA,T)
   call MergeSort(A(NA+1),NB,T)
 
   if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      call Merge(T,NA,A(NA+1),NB,A,N)
   endif
   return
 
end subroutine MergeSort
 
program TestMergeSort
 
   integer, parameter :: N = 14
   integer, dimension(N) :: A = (/ 1, 5, 45, 2, 7, 3, 9, 4, 6, 10, 11, 22, 13, 14 /)
   integer, dimension ((N+1)/2) :: T
   integer :: ip

   call MergeSort(A,N,T)

   write(*,*) 'Sorted array:'
   do ip = 1,N
     write(*,'(1x,I5,I8)') ip, A(ip)
   enddo
 
end program

