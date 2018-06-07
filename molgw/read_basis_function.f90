!=========================================================================
subroutine read_basis_function(unitfile,bf)
  use m_basis_set
 implicit none

 integer,intent(in)               :: unitfile
 type(basis_function),intent(out) :: bf
!=====
!=====

 read(unitfile)  bf%shell_index
 read(unitfile)  bf%am
 read(unitfile)  bf%amc
 read(unitfile)  bf%nx
 read(unitfile)  bf%ny
 read(unitfile)  bf%nz
 read(unitfile)  bf%iatom
 read(unitfile)  bf%x0(:)
 read(unitfile)  bf%ngaussian
 allocate(bf%g(bf%ngaussian))
 read(unitfile)  bf%g(:)
 allocate(bf%coeff(bf%ngaussian))
 read(unitfile)  bf%coeff(:)

end subroutine read_basis_function