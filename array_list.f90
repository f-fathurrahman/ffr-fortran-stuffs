MODULE array_list

  TYPE container
    class(*), allocatable :: item
    class(*), allocatable :: items(:)
  END TYPE

  INTERFACE add_item
    MODULE PROCEDURE add_item_scalar
    MODULE PROCEDURE add_item_array
  END INTERFACE

CONTAINS 

  SUBROUTINE add_item_scalar(a, e)
    TYPE(container), ALLOCATABLE, INTENT(inout) :: a(:)
    class(*), INTENT(in) :: e
    TYPE(container),ALLOCATABLE :: tmp(:)

      if (.not.allocated(a)) then
        allocate(a(1))
        allocate(a(1)%item, source = e)
      else
        call move_alloc(a,tmp)
        allocate(a(size(tmp)+1))
        a(1:size(tmp)) = tmp
        allocate(a(size(tmp)+1)%item, source = e)
      end if
   end SUBROUTINE

  subroutine add_item_array(a, e)
    type(container),allocatable,intent(inout) :: a(:)
    class(*),intent(in) :: e(:)
    type(container),allocatable :: tmp(:)

      if (.not.allocated(a)) then
        allocate(a(1))
        allocate(a(1)%items(size(e)), source = e)
      else
        call move_alloc(a,tmp)
        allocate(a(size(tmp)+1))
        a(1:size(tmp)) = tmp
        allocate(a(size(tmp)+1)%items(size(e)), source = e)
      end if
   end SUBROUTINE

end module



  use array_list

  type(container), allocatable :: a_list(:)

  type newtype
  end TYPE

  integer i

  call add_item(a_list, 1)

  call add_item(a_list, 5.5)

  call add_item(a_list, (4., 5.))

  call add_item(a_list, newtype())

  call add_item(a_list, [1, 2, 3, 4, 5])

  do i = 1, size(a_list)
    call print(a_list(i))
  end DO

contains

  subroutine print(c)
    type(container), intent(in) :: c

    if (allocated(c%item)) then
      select type (x=>c%item)
        type is (integer)
          print *, x
        type is (real)
          print *, x
        type is (complex)
          print *, x
        type is (newtype)
          print *, "is newtype"
        class default
          print *, "is unknown type"
      end select
    else if (allocated(c%items)) then
      select type (x=>c%items)
        type is (integer)
          print *, x
        type is (real)
          print *, x
        type is (complex)
          print *, x
        type is (newtype)
          print *, "is array of newtype"
        class default
          print *, "is array of unknown type"
      end select
    else
      write(*,*) "Error, empty container"
    end if
  END SUBROUTINE

END 

