program test
 use iso_c_binding
 implicit none
 interface
  subroutine testc(pa) bind(c)
   use iso_c_binding
   type(c_ptr):: pa
  end subroutine testc
 end interface
 
 type(c_ptr) :: pa
 real(c_double),pointer::fpa(:)
 
 call testc(pa)
 call c_f_pointer(pa, fpa, [5])
 print*, fpa(1)
 
end program test

