#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program abinit

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_parameters
 use m_ab7_invars
 use m_build_info
 use m_cppopts_dumper
 use m_optim_dumper
 use m_profiling_abi
 use m_results_out
 use m_xmpi
 use m_xomp
 use m_xpapi
 use m_errors
 use m_argparse
 use m_nctk
#if defined HAVE_MPI2
 use mpi
#endif

 use defs_time,     only : time_set_papiopt
 use m_time ,       only : asctime, sec2str
 use m_fstrings,    only : sjoin, strcat, itoa, yesno, ljust
 use m_io_tools,    only : open_file, flush_unit, delete_file, num_opened_units, show_units
 use m_specialmsg,  only : specialmsg_getcount
 use m_exit,        only : get_timelimit_string
 use m_atomdata,    only : znucl2symbol
 use m_libpaw_tools,only : libpaw_spmsg_getcount
 use m_pawxmlps,    only : paw_setup, paw_setup_free, npsp_pawxml,ipsp2xml
 use m_mpinfo,      only : destroy_mpi_enreg, clnmpi_img, clnmpi_grid, clnmpi_atom, clnmpi_pert
#ifdef HAVE_GPU_CUDA
 use m_initcuda,     only: setdevice_cuda,unsetdevice_cuda
#endif
#if defined HAVE_BIGDFT
 use BigDFT_API,    only : bigdft_init_errors,bigdft_init_timing_categories
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abinit'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_57_iovars
 use interfaces_95_drive
!End of the abilint section

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments -----------------------------------
!Local variables-------------------------------
!
!===============================================================================
!  abinit_version designate overall code version
!  mpw=maximum number of planewaves in basis sphere
!  unit numbers (ab_in,ab_out,std_out,tmp_unit) have been defined in defs_basis.f .
!  The array filnam is used for the name of input and output files,
!  and roots for generic input, output or temporary files.
!  Pseudopotential file names are set in iofn2, and are contained in pspheads.
!  The name filstat will be needed beyond gstate to check
!  the appearance of the "exit" flag, to make a hasty exit, as well as
!  in order to output the status of the computation.
!==============================================================================
! Declarations
! Define "level of the routine", for debugging purposes
 integer,parameter :: level=1
 integer :: choice,dmatpuflag,ierr,iexit,ii,iounit,dtsetsId,ios
 integer :: lenstr,me,print_mem_report
 integer :: mu,natom,ncomment,ncomment_paw,ndtset
 integer :: ndtset_alloc,nexit,nexit_paw,nfft,nkpt,npsp
 integer :: nsppol,nwarning,nwarning_paw,papiopt,prtvol,timopt,use_gpu_cuda
 integer,allocatable :: nband(:),npwtot(:)
 real(dp) :: etotal
 real(dp) :: tcpui,twalli
 real(dp) :: strten(6),tsec(2)
 real(dp),allocatable :: fred(:,:),xred(:,:)
 character(len=24) :: codename
 character(len=24) :: start_datetime
 character(len=5000) :: message
 character(len=strlen) :: string
 character(len=fnlen) :: filstat
 character(len=fnlen) :: filnam(5)
 type(args_t) :: args
 type(dataset_type),pointer  :: dtsets(:)
 type(MPI_type),pointer :: mpi_enregs(:)
 type(pspheader_type),pointer :: pspheads(:)
 type(results_out_type),allocatable,target :: results_out(:)
 type(results_out_type),pointer :: results_out_all(:)
 type(ab_dimensions) :: mxvals
 logical :: test_img,test_exit,use_results_all,xml_output=.false.
 integer :: values(8)
 character(len=5) :: strzone
 character(len=8) :: strdat
 character(len=10) :: strtime
 character(len=13) :: warn_fmt
#ifdef HAVE_GPU_CUDA
 integer :: gpu_devices(5)
#endif
#if defined HAVE_MPI
 real(dp) :: tsec_s(2)
#endif


!******************************************************************

!0) Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 !call xlf_set_sighandler()

!------------------------------------------------------------------------------

!1) Eventually initialize MPI
!Pay attention: it may be initialzed again in finddistrproc

 call xmpi_init()
 me=xmpi_comm_rank(xmpi_world)

 ! Parse command line arguments.
 args = args_parser(); if (args%exit /= 0) goto 100

!Initialize memory profiling if it is activated
!if a full memocc.prc report is desired, set the argument of abimem_init to "2" instead of "0"
!note that memocc.prc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(args%abimem_level)
!call abimem_init(2)
#endif

!------------------------------------------------------------------------------

!2) Initialize overall timing of run:
 call xpapi_init()
 call xpapi_show_info(unit=std_out,mode_paral="COLL")

 start_datetime = asctime()

 call timein(tcpui,twalli)

 call timab(1,0,tsec)

!Start to accumulate time for the entire run. The end of accumulation is in timana.f
 call timab(1,1,tsec)

!------------------------------------------------------------------------------

!3) Print greeting for interactive user,
!read names of files (input, output, rootinput, rootoutput, roottemporaries),
!create the name of the status file, initialize the status subroutine.

 call timab(41,3,tsec)

 call iofn1(filnam,filstat,xmpi_world)

!------------------------------------------------------------------------------

!4) Open output file and print herald at top of output and log files

 if (me==0) then
#ifdef FC_NAG
   open(unit=ab_out,file=filnam(2),form='formatted',status='new', action="write", recl=ABI_RECL, iomsg=message, iostat=ios)
#else
   open(unit=ab_out,file=filnam(2),form='formatted',status='new', action="write", iomsg=message, iostat=ios)
#endif
   
   ABI_CHECK(ios == 0, message)
   
   rewind (unit=ab_out)
   
   codename='ABINIT'//repeat(' ',18)

   call herald(codename,abinit_version,ab_out)
   call herald(codename,abinit_version,std_out)

   call dump_config(std_out)
   call dump_optim(std_out)
   call dump_cpp_options(std_out)
!  Write names of files
   write(message, '(a,a,a,a,a,a,a,a,a,a,a,a)' )&
&   '- input  file    -> ',trim(filnam(1)),ch10,&
&   '- output file    -> ',trim(filnam(2)),ch10,&
&   '- root for input  files -> ',trim(filnam(3)),ch10,&
&   '- root for output files -> ',trim(filnam(4)),ch10

      WRITE(*,*) 'ab_out = ', ab_out
      WRITE(*,*) 'std_out = ', std_out
      WRITE(*,*) 'message = ', trim(message)


   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 ! Test if the netcdf library supports MPI-IO
 call nctk_test_mpiio()






  IF( me == 0 ) CLOSE(unit=ab_out)

100 CALL xmpi_end()

end program 


