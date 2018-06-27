program segfault1
  implicit none
  call foo(4096000)
contains
  subroutine foo(n)
    integer, intent(in) :: n
    integer, dimension(n) :: a, b
    a = 1
    b = 2
    a = a + b
    print *, sum(a)
  end subroutine foo
end program segfault1
