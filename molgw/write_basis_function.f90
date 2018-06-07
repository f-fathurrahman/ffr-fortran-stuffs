!=========================================================================
subroutine write_basis_function(unitfile,bf)
  use m_basis_set
 implicit none

 integer,intent(in)              :: unitfile
 type(basis_function),intent(in) :: bf
!=====
!=====

 write(unitfile)  bf%shell_index  
 write(unitfile)  bf%am           
 write(unitfile)  bf%amc          
 write(unitfile)  bf%nx
 write(unitfile)  bf%ny
 write(unitfile)  bf%nz
 write(unitfile)  bf%iatom        
 write(unitfile)  bf%x0(:)        
 write(unitfile)  bf%ngaussian    
 write(unitfile)  bf%g(:)         
 write(unitfile)  bf%coeff(:)     


end subroutine write_basis_function