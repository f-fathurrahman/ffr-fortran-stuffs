! from: https://sites.google.com/site/rupakmukherjee01/program-archive/fortran-95/basic-programs/fftw-2-d-serial-real-to-complex-to-real

Program rupak

implicit none

integer ( kind = 4 ), parameter :: Nx = 5
integer ( kind = 4 ), parameter :: Ny = 3
integer ( kind = 4 ), parameter :: Nh = ( Nx / 2 ) + 1

include "fftw3.f"

integer ( kind = 4 ) i,j

complex ( kind = 8 ) ux(Nh,Ny)

real ( kind = 8 ) ukx(Nx,Ny),ukx_dum(Nx,Ny)

integer ( kind = 8 ) plan_forward,plan_backward

  integer :: seed

  seed = 123456789

!integer,parameter :: seed = 99999999

! call srand(seed)

 

do i = 1,Nx

  do j = 1,Ny

    ukx(i,j) = ran(seed)!rand()!(1.0d0,0.0d0)*rand() + (0.0d0,1.0d0)*rand()

     seed=ukx(i,j)*seed

     write(*, '(2x,i4,2x,i4,2x,2g14.6)') i,j,ukx(i,j), "ukx"

  enddo

enddo

  call dfftw_plan_dft_r2c_2d_ (plan_forward, Nx, Ny, ukx, ux, FFTW_ESTIMATE)

  call dfftw_execute_ (plan_forward)

do i = 1,Nx/2+1

  do j = 1,Ny

    ux(i,j) = ux(i,j)

    !write(*, '(2x,i4,2x,i4,2x,2g14.6)') i,j,ux(i,j)

  enddo

enddo

  call dfftw_plan_dft_c2r_2d_ (plan_backward, Nx, Ny, ux, ukx_dum, FFTW_ESTIMATE)

  call dfftw_execute_ (plan_backward)

do i = 1,Nx

  do j = 1,Ny

    write(*, '(2x,i4,2x,i4,2x,2g14.6)') i,j,ukx_dum(i,j)/(float(Nx)*float(Ny)), "ukx_dum"

  enddo

enddo

end