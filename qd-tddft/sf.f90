!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROGRAM STRENGTH_FUNCTION
! =========================
!
! This program calculates the sine Fourier transform of the electronic dipole
! written by "td" in the "dipole" file. The result is written in "spectrum"
! file.
!
! The sine Fourier transform is previously filtered in time, by multiplying
! the original function by a third-order polynomial that is one at time zero,
! and zero at the final time.
!
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine strength_function
  implicit none

  integer :: i, j, unit, nsteps
  real(8) :: x, dt, e, t, relt, sf(2), filter
  real(8), allocatable :: signal(:, :)

  character(len=8), parameter :: infilename  = "dipole", &
                                 outfilename = "spectrum"
  integer, parameter :: nesteps = 2001
  real(8), parameter :: dw = 0.001_8, &
                        emax = 0.001_8*(nesteps-1)

  unit = 11
  open(unit=unit, file=trim(infilename), status="old", action="read", iostat = i)
  if(i.ne.0) then
     write(0, '(a)') "Error opening file '"//trim(infilename)//"'."; stop
  endif

  nsteps = 0
  do
    read(unit, *, iostat = i) x, x, x
    if(i.ne.0) exit
    nsteps = nsteps + 1
  enddo

  if(nsteps.eq.0) then
     write(0,'(a)') "File '"//trim(infilename)//"' seems to be empty or unreadable."; stop
  endif

  rewind(unit)
  read(unit, *) x
  read(unit, *) dt
  rewind(unit)

  allocate(signal(nsteps, 2))

  do i = 1, nsteps
     read(unit, *) x, signal(i, 1), signal(i, 2)
  enddo

  close(unit)
  open(unit=unit, file=trim(outfilename), status="replace", action="write", iostat = i)

  do i = 1, nesteps
     e = (i-1)*dw
     sf = 0.0_8
     do j = 1, nsteps
        t = (j-1)*dt; relt = (real(j-1, 8)/real(nsteps-1, 8))
        filter = 1.0_8 - 3.0_8*relt**2 + 2.0_8*relt**3
        sf(1) = sf(1) + sin(e*t)*signal(j, 1)*filter
        sf(2) = sf(2) + sin(e*t)*signal(j, 2)*filter
     enddo
     write(unit, '(3e18.8)') e, sf(1), sf(2)
  enddo

  close(unit = unit)
end subroutine strength_function
