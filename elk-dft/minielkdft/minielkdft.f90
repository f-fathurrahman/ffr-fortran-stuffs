
! Copyright (C) 2002-2011 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! main routine for the Elk code
program elk
use modmain
use modmpi
use modomp
use modvars
implicit none
! local variables
logical exist
integer itask
! initialise MPI execution environment
call mpi_init(ierror)
! duplicate mpi_comm_world
call mpi_comm_dup(mpi_comm_world,mpicom,ierror)
! determine the number of MPI processes
call mpi_comm_size(mpicom,np_mpi,ierror)
! determine the local MPI process number
call mpi_comm_rank(mpicom,lp_mpi,ierror)
! determine if the local process is the master
if (lp_mpi.eq.0) then
  mp_mpi=.true.
  write(*,*)
  write(*,'("Elk code version ",I1.1,".",I1.1,".",I2.2," started")') version
else
  mp_mpi=.false.
end if

! read input files
call readinput

! initialise OpenMP variables
call omp_init

! initialise the MKL library
call mkl_init

! initialise the OpenBLAS library
call oblas_init

! initialise the BLIS library
call blis_init

if (mp_mpi) then
  write(*,*)
  write(*,'("Number of MPI processes : ",I6)') np_mpi
  write(*,'("Number of OpenMP threads per MPI process : ",I4)') maxthd
  write(*,'("Total number of threads : ",I6)') np_mpi*maxthd
  write(*,'("Maximum OpenMP nesting level : ",I4)') maxlvl
end if

! delete the VARIABLES.OUT file
call delvars

! write version number to VARIABLES.OUT
call writevars('version',nv=3,iva=version)

! check if Elk is already running in this directory
if (mp_mpi) then
  inquire(file='RUNNING',exist=exist)
  if (exist) then
    write(*,*)
    write(*,'("Info(elk): several copies of Elk may be running in this path")')
    write(*,'("(this could be intentional, or result from a previous crash,")')
    write(*,'(" or arise from an incorrect MPI compilation)")')
  else
    open(50,file='RUNNING')
    close(50)
  end if
end if

! perform the tasks
do itask=1,ntasks
  task=tasks(itask)
  if (mp_mpi) then
    write(*,*)
    write(*,'("Info(elk): current task : ",I6)') task
  end if

! synchronise MPI processes
  call mpi_barrier(mpicom,ierror)

! check if task can be run with MPI
  if (lp_mpi.gt.0) then
    if (any(task.eq.[0,1,2,3,5,15,16,28,29,61,62,63,110,120,135,136,162,170, &
     180,185,188,200,201,205,240,241,270,300,320,330,331,350,351,360,371,372, &
     373,440,460,461,600,610,620,630,640,700,701])) then
      continue
    else
      write(*,'("Info(elk): MPI process ",I6," idle for task ",I6)') lp_mpi,task
      cycle
    end if
  end if

! write task to VARIABLES.OUT
  call writevars('task',iv=task)
  select case(task)
  case(0,1)
    call gndstate
!  case(2,3)
!    call geomopt
!  case(5)
!    call hartfock
!  case(10)
!    call writedos
!  case(14)
!    call writesf
!  case(15,16)
!    call writelsj
!  case(20,21,22,23)
!    call bandstr
!  case(25)
!    call effmass
!  case(28,29)
!    call mae
!  case(31,32,33)
!    call rhoplot
!  case(41,42,43)
!    call potplot
!  case(51,52,53)
!    call elfplot
!  case(61,62,63,162)
!    call wfplot
!  case(65)
!    call wfcrplot
!  case(71,72,73,81,82,83,141,142,143,151,152,153)
!    call vecplot
!  case(91,92,93)
!    call dbxcplot
!  case(100,101)
!    call fermisurf
!  case(102)
!    call fermisurfbxsf
!  case(105)
!    call nesting
!  case(110)
!    call mossbauer
!  case(115)
!    call writeefg
!  case(120)
!    call writepmat
!  case(121)
!    call dielectric
!  case(122)
!    call moke
!  case(125)
!    call nonlinopt
!  case(130)
!    call writeexpmat
!  case(135)
!    call writewfpw
!  case(140)
!    call elnes
!  case(150)
!    call writeevsp
!  case(160)
!    call torque
!  case(170)
!    call writeemd
!  case(171,172,173)
!    call emdplot
!  case(180)
!    call writeepsinv
!  case(185)
!    call writehmlbse
!  case(186)
!    call writeevbse
!  case(187)
!    call dielectric_bse
!  case(190)
!    call geomplot
!  case(195)
!    call sfacrho
!  case(196)
!    call sfacmag
!  case(200,201,202)
!    call phononsc
!  case(205)
!    call phonon
!  case(210)
!    call phdos
!  case(220)
!    call phdisp
!  case(230)
!    call writephn
!  case(240,241)
!    call ephcouple
!  case(245)
!    call phlwidth
!  case(250)
!    call alpha2f
!  case(260)
!    call eliashberg
!  case(270)
!**************
!   ! ?????
!  case(300)
!    call rdmft
!  case(320)
!    call tddftlr
!  case(330,331)
!    call tddftsplr
!  case(341,342,343)
!    call wxcplot
!  case(350,351,352)
!    call spiralsc
!  case(360)
!    call ssfsmjx
!  case(371,372,373)
!    call curdenplot
!  case(400)
!    call writetmdu
!  case(430)
!    call writestrain
!  case(440)
!    call writestress
!  case(450)
!    call genafieldt
!  case(460,461)
!    call tddft
!  case(480,481)
!    call dielectric_tdrt
!  case(500)
!    call testcheck
!  case(550)
!    call writew90
!  case(600)
!    call gwsefm
!  case(610)
!    call gwspecf
!  case(620)
!    call gwbandstr
!  case(630)
!    call gwscrho
!  case(640)
!    call gwdmat
!  case(700,701)
!    call gndstulr
!  case(731,732,733)
!    call rhouplot
!  case(741,742,743)
!    call potuplot
!  case(771,772,773)
!    call maguplot
  case default
    write(*,*)
    write(*,'("Error(elk): task not defined : ",I8)') task
    write(*,*)
    stop
  end select
! reset the OpenMP thread variables
  call omp_reset
! close all opened files
  call closefiles
end do
if (mp_mpi) then
  open(50,file='RUNNING')
  close(50,status='DELETE')
  write(*,*)
  write(*,'("Elk code stopped")')
end if
! terminate MPI execution environment
call mpi_finalize(ierror)
end program

