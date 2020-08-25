program version_options                    
  use iso_fortran_env, only : compiler_version, compiler_options
  implicit none                                                                         
  write(*,*)"Compiler: ", compiler_version()
  write(*,*)"Options: ", compiler_options()
end program version_options