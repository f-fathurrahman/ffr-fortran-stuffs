subroutine header()
#ifdef FORTRAN2008
 use,intrinsic :: iso_fortran_env, only: compiler_version,compiler_options
#endif
 use m_definitions
 use m_mpi
 use m_warning,only: issue_warning
 use m_tools,only: orbital_momentum_name
 use m_libint_tools,only: libint_init
 implicit none

#ifdef _OPENMP
 integer,external  :: OMP_get_max_threads
 character(len=64) :: msg
#endif
 character(len=40)   :: git_sha
 integer             :: values(8) 
 character(len=1024) :: chartmp
 integer             :: nchar,kchar,lchar
!=====
! variables used to call C
 integer(C_INT)      :: ammax
 logical(C_BOOL)     :: has_onebody,has_gradient
!=====

! Here we call the fortran code that was generated by the python script
! Any new variable should be added through the python script
#include "INCLUDED/git_sha.f90"
 
!=====

 write(stdout,'(1x,70("="))') 
 write(stdout,'(/,/,12x,a,/)') 'Welcome to the fascinating world of MOLGW'
 write(stdout,'(24x,a)')       'version 1.F'
 write(stdout,'(/,/,1x,70("="))') 

 write(stdout,'(/,a,a,/)') ' MOLGW commit git SHA: ',git_sha
#ifdef FORTRAN2008
 write(stdout,'(1x,a,a)')    'compiled with ',compiler_version()
 write(stdout,'(1x,a)')      'with options: '
 chartmp = compiler_options()
 nchar = LEN(TRIM(chartmp))
 kchar = 1
 lchar = 0
 do 
   lchar = SCAN(chartmp(kchar:nchar),' ')
   if( lchar == 0 ) exit
   write(stdout,'(6x,a,a)') 'FCOPT  ',chartmp(kchar:kchar+lchar-1)
   kchar = kchar + lchar
 enddo
 write(stdout,*)
#endif


 call date_and_time(VALUES=values)

 write(stdout,'(a,i2.2,a,i2.2,a,i4.4)') ' Today is ',values(2),'/',values(3),'/',values(1)
 write(stdout,'(a,i2.2,a,i2.2)')        ' It is now ',values(5),':',values(6)
 select case(values(5))
 case(03,04,05,06,07)
   write(stdout,*) 'And it is too early to work. Go back to sleep'
 case(22,23,00,01,02)
   write(stdout,*) 'And it is too late to work. Go to bed and have a sleep'
 case(12,13)
   write(stdout,*) 'Go and get some good food'
 case(17)
   write(stdout,*) 'Dont forget to go and get the kids'
 case default
   write(stdout,*) 'And it is perfect time to work'
 end select


 write(stdout,'(/,1x,a)') 'Linking options:'
#ifdef HAVE_LIBXC
!#ifndef LIBXC_SVN
! call xc_f90_version(values(1),values(2))
! write(chartmp,'(i2,a,i2)') values(1),'.',values(2)
!#else
! call xc_f90_version(values(1),values(2),values(3))
! write(chartmp,'(i2,a,i2,a,i2)') values(1),'.',values(2),'.',values(3)
!#endif
! write(stdout,*) 'LIBXC version '//TRIM(chartmp)
 write(stdout,*) 'Running with LIBXC'
#endif
#ifdef _OPENMP
 write(msg,'(i6)') OMP_get_max_threads()
 msg='OPENMP option is activated with threads number'//msg
 call issue_warning(msg)
#endif
#if defined HAVE_MPI && defined HAVE_SCALAPACK
 write(stdout,*) 'Running with MPI'
 write(stdout,*) 'Running with SCALAPACK'
#endif
#if defined(HAVE_MPI) && !defined(HAVE_SCALAPACK)
 call die('Code compiled with SCALAPACK, but without MPI. This is not permitted')
#endif
#if !defined(HAVE_MPI) && defined(HAVE_SCALAPACK)
 call die('Code compiled with MPI, but without SCALAPACK. This is not permitted')
#endif

 ! LIBINT details
 call libint_init(ammax,has_onebody,has_gradient)
 write(stdout,'(1x,a)')        'Running with LIBINT (to calculate the Coulomb integrals)'
 write(stdout,'(6x,a,i5,3x,a)') 'max angular momentum handled by your LIBINT compilation: ', &
                                ammax,orbital_momentum_name(ammax)
 call set_molgw_lmax(ammax)

#ifdef HAVE_LIBINT_ONEBODY
 if( .NOT. has_onebody ) &
   call die('MOLGW compiled with LIBINT one-body terms, however the LIBINT compilation does not calculate the one-body terms')
 if( .NOT. has_gradient ) &
   call die('LIBINT compilation does not have the first derivative')
 write(stdout,'(1x,a)') 'Using LIBINT for the one-body parts of the Hamiltonian as well'
#endif
 write(stdout,*)
 write(stdout,*)


end subroutine header
