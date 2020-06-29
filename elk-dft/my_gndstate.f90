subroutine my_gndstate()
  ! initialise global variables
  call init0()
  call init1()
  
  ! initialise the density and magnetisation from atomic data
  call rhoinit()
  call maginit()
  
  ! generate the Kohn-Sham potential and magnetic field
  call potks(.true.)
end subroutine