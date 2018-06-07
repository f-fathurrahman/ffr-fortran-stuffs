!=========================================================================
subroutine dump_out_matrix(print_matrix,title,n,nspin,matrix)
 use m_definitions
 use m_mpi
 implicit none
 logical,intent(in)          :: print_matrix       
 character(len=*),intent(in) :: title
 integer,intent(in)          :: n,nspin
 real(dp),intent(in)         :: matrix(n,n,nspin)
!=====
 integer,parameter :: MAXSIZE=25
!=====
 integer :: i,ispin
!=====

 if( .NOT. print_matrix ) return

 write(stdout,'(/,1x,a)') TRIM(title)

 do ispin=1,nspin
   if(nspin==2) then
     write(stdout,'(a,i1)') ' spin polarization # ',ispin
   endif
   do i=1,MIN(n,MAXSIZE)
     write(stdout,'(1x,i3,100(1x,f12.5))') i,matrix(i,1:MIN(n,MAXSIZE),ispin)
   enddo
   write(stdout,*)
 enddo
 write(stdout,*)

end subroutine dump_out_matrix