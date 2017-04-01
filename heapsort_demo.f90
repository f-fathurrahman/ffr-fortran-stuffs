program Heapsort_Demo
  implicit none
 
  integer, parameter :: N = 5
  real :: array(N)
  INTEGER :: ip
 
  call random_seed
  call random_number(array)
  write(*,*) "Unsorted array:-"
  DO ip = 1, N
    write(*,'(1x,I5,F18.10)') ip, array(ip)
  ENDDO
 
  write(*,*)
  call heapsort(array)
  write(*,*) "Sorted array:-"
  DO ip = 1, N
    write(*,'(1x,I5,F18.10)') ip, array(ip)
  ENDDO

CONTAINS 
 
SUBROUTINE heapsort(a)
 
   real, intent(in out) :: a(0:)
   integer :: start, n, bottom
   real :: temp
 
   n = size(a)
   WRITE(*,*) 'N = ', N
   do start = (n - 2) / 2, 0, -1
     WRITE(*,*) 'start = ', start
     call siftdown(a, start, n);
   end DO
 
   WRITE(*,*)
   do bottom = n - 1, 1, -1
     WRITE(*,*) 'bottom = ', bottom
     temp = a(0)
     a(0) = a(bottom)
     a(bottom) = temp;
     call siftdown(a, 0, bottom)
   end do
 
end subroutine heapsort
 
subroutine siftdown(a, start, bottom)
 
  real, intent(in out) :: a(0:)
  integer, intent(in) :: start, bottom
  integer :: child, root
  real :: temp
 
  root = start
  do while(root*2 + 1 < bottom)
    child = root * 2 + 1
 
    if (child + 1 < bottom) then
      if (a(child) < a(child+1)) child = child + 1
    end if
 
    if (a(root) < a(child)) then
      temp = a(child)
      a(child) = a (root)
      a(root) = temp
      root = child
    else
      return
    end if  
  end do      
 
end subroutine siftdown
 
end program Heapsort_Demo
