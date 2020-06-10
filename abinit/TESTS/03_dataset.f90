#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

PROGRAM abinit
  
  USE m_argparse, ONLY : args_t
  USE m_argparse, ONLY : args_parser

  USE defs_abitypes, ONLY : MPI_type, dataset_type, ab_dimensions
  USE defs_datatypes, ONLY : pspheader_type
  USE defs_basis, ONLY : dp, fnlen, ab_out, std_out, ch10, strlen
  USE defs_basis, ONLY : abi_io_redirect

  USE defs_time, ONLY : time_set_papiopt
  
  USE m_xmpi, ONLY : xmpi_world
  USE m_xmpi, ONLY : xmpi_init, xmpi_comm_rank, xmpi_end

  USE m_build_info, ONLY : abinit_version
  USE m_build_info, ONLY : dump_config
  USE m_optim_dumper, ONLY : dump_optim
  USE m_cppopts_dumper, ONLY : dump_cpp_options
  USE m_nctk, ONLY : nctk_test_mpiio

  USE m_results_out

  USE m_xpapi, ONLY : xpapi_init, xpapi_show_info
  USE m_time, ONLY : asctime

  USE m_ab7_invars
  USE m_errors

  USE interfaces_14_hidewrite

  IMPLICIT NONE 

  INTEGER :: me
  TYPE(args_t) :: args
  CHARACTER(len=fnlen) :: filnam(5)
  CHARACTER(len=fnlen) :: filstat
  CHARACTER(len=5000) :: message
  REAL(DP) :: tsec(2), tcpui, twalli
  CHARACTER(len=24) :: start_datetime
  INTEGER :: ndtset
  CHARACTER(len=strlen) :: string
  INTEGER :: lenstr

  TYPE(MPI_type), POINTER :: mpi_enregs(:)

  INTEGER :: dmatpuflag

  TYPE(dataset_type), POINTER :: dtsets(:)
  INTEGER :: dtsetsId, ndtset_alloc

  INTEGER :: npsp
  INTEGER :: papiopt
  INTEGER :: timopt

  TYPE(ab_dimensions) :: mxvals

  TYPE(pspheader_type), POINTER :: pspheads(:)

  INTEGER :: ios
  CHARACTER(len=24) :: codename

  INTEGER :: iexit
  INTEGER :: level=1 ! "level of the routine", for debugging purpose

  TYPE(results_out_type), ALLOCATABLE, TARGET :: results_out(:)
  TYPE(results_out_type), POINTER :: results_out_all(:)

  INTEGER :: choice
  LOGICAL :: test_img, use_results_all

  INTEGER :: iounit
  INTEGER :: values(8)
  LOGICAL :: xml_output = .false.



  CALL abi_io_redirect(new_io_comm=xmpi_world)

  CALL xmpi_init()
  me = xmpi_comm_rank(xmpi_world)

  ! Parse command line arguments.
  args = args_parser()
  IF( args%EXIT /= 0 ) GOTO 100

  CALL step_02()

  CALL step_03()

  CALL step_04()

  CALL step_05()

  CALL step_06_11()

  CALL step_12()

  CALL test_dataset()

  IF( me == 0 ) CLOSE(unit=ab_out)

100 CALL xmpi_end()


CONTAINS 


  SUBROUTINE step_02()

    WRITE(*,*) '-------'
    WRITE(*,*) 'Step 02'
    WRITE(*,*) '-------'

    CALL xpapi_init()
    CALL xpapi_show_info(unit=std_out,mode_paral="COLL")

    start_datetime = asctime()

    CALL timein(tcpui,twalli)

    CALL timab(1,0,tsec)

    !Start to accumulate time for the entire run. The end of accumulation is in timana.f
    CALL timab(1,1,tsec)

  END SUBROUTINE 


  SUBROUTINE step_03()

    WRITE(*,*) '-------'
    WRITE(*,*) 'Step 03'
    WRITE(*,*) '-------'

    CALL timab( 41, 3, tsec )
    ! Greetings to interactive user
    CALL iofn1( filnam, filstat, xmpi_world )
  END SUBROUTINE 


  SUBROUTINE step_04()

    WRITE(*,*) '-------'
    WRITE(*,*) 'Step 04'
    WRITE(*,*) '-------'

    IF( me==0 ) THEN 

      OPEN( unit=ab_out, file=filnam(2),form='formatted',status='new', &
            action='write', iomsg=message, iostat=ios)

      ABI_CHECK(ios == 0, message)

      REWIND( unit=ab_out )

      codename = 'ABINIT'//repeat(' ',18)
      
      CALL herald( codename, abinit_version, ab_out )
      CALL herald( codename, abinit_version, std_out )
      
      CALL dump_config( std_out )
      CALL dump_optim( std_out )
      CALL dump_cpp_options( std_out )

      WRITE(*,*) '-------- After dump_cpp_options -------'
      
      !
      ! Write names of files
      !
      WRITE(message, '(a,a,a,a,a,a,a,a,a,a,a,a)' ) &
          &   '- input  file    -> ',trim(filnam(1)),ch10,&
          &   '- output file    -> ',trim(filnam(2)),ch10,&
          &   '- root for input  files -> ',trim(filnam(3)),ch10,&
          &   '- root for output files -> ',trim(filnam(4)),ch10
      WRITE(*,*) trim(filnam(1))
      WRITE(*,*) trim(filnam(2))
      WRITE(*,*) trim(filnam(3))
      WRITE(*,*) trim(filnam(4))
      WRITE(*,*) 'ab_out = ', ab_out
      WRITE(*,*) 'std_out = ', std_out
      WRITE(*,*) 'message = ', trim(message)

      CALL wrtout( ab_out, message, 'COLL' )
      CALL wrtout( std_out, message, 'COLL' )
    ENDIF 

    ! Test if the netcdf library supports MPI-IO
    CALL nctk_test_mpiio()

  END SUBROUTINE 



  !5) Read the file, stringify it and return the number of datasets.
  SUBROUTINE step_05()

    WRITE(*,*) '-------'
    WRITE(*,*) 'Step 05'
    WRITE(*,*) '-------'

    CALL parsefile( filnam(1), lenstr, ndtset, string, xmpi_world )
  END SUBROUTINE 

  
  SUBROUTINE step_06_11()
    
    WRITE(*,*) '----------'
    WRITE(*,*) 'Step 06-11'
    WRITE(*,*) '----------'

    !6~11) Call the parser from the parser module.
    CALL ab7_invars_set_flags(.true., .true., status_file = filstat, timab_tsec = tsec)
    CALL ab7_invars_load(dtsetsId, string, lenstr, ndtset, .true., .true.)
    CALL timab(44,1,tsec)
    CALL ab7_invars_get_abinit_vars(dtsetsId, dtsets, pspheads,mxvals, papiopt, timopt, dmatpuflag)
    
    ndtset_alloc = size(dtsets) - 1
    npsp = size(pspheads)

    !Enable PAPI timers
    CALL time_set_papiopt(papiopt)

    ABI_DATATYPE_ALLOCATE(mpi_enregs,(0:max(1,ndtset)))
    CALL mpi_setup(dtsets,filnam,lenstr,mpi_enregs,ndtset,ndtset_alloc,string)

    CALL memory_eval(dtsets, ab_out, mpi_enregs, ndtset, ndtset_alloc, npsp, pspheads)
  END SUBROUTINE 


  SUBROUTINE step_12()

    INTEGER :: ii
    CHARACTER(len=5) :: strzone
    CHARACTER(len=8) :: strdat
    CHARACTER(len=10) :: strtime

    WRITE(*,*) '-------'
    WRITE(*,*) 'Step 12'
    WRITE(*,*) '-------'
    
    CALL status(0, filstat, iexit, level, 'call outvars(1)')

    ABI_DATATYPE_ALLOCATE( results_out, (0:ndtset_alloc) )

    CALL init_results_out( dtsets, 1, 1, mpi_enregs, &
                           mxvals%natom, &
                           mxvals%mband_upper, &
                           mxvals%nkpt, &
                           npsp, &
                           mxvals%nsppol, &
                           mxvals%ntypat, &
                           results_out )

    ! Gather contributions to results_out from images of the cell, if needed
    test_img = ( mxvals%nimage /= 1 .and. maxval(dtsets(:)%npimage) > 1 )
    use_results_all=.false.
    
    IF( test_img ) THEN 
      
      use_results_all = ( me == 0 )
      
      IF( use_results_all )  THEN 
         ABI_DATATYPE_ALLOCATE(results_out_all,(0:ndtset_alloc))
      ENDIF

      CALL gather_results_out( dtsets, mpi_enregs, results_out, results_out_all, use_results_all, & 
                               allgather=.false., master=0 )

    ELSE
      
      WRITE(*,*) 'ffr: Not test_img'
      results_out_all => results_out
    
    ENDIF 

    
    IF( me == 0 ) THEN 
      
      !  Echo input to output file on unit ab_out, and to log file on unit 06 :
      choice=1
      
      DO ii = 1,2
        IF( ii == 1 ) iounit = ab_out
        IF( ii == 2 ) iounit = std_out
        CALL outvars( choice, dmatpuflag, dtsets, trim(filnam(4)), &
                      iounit, mxvals, ndtset, ndtset_alloc, npsp, results_out_all, timopt )
      ENDDO

      IF( dtsets(1)%prtxml == 1) THEN 
        CALL outxml_open(trim(filnam(4)))
        CALL date_and_time( strdat, strtime, strzone, values )
        xml_output = .true.
      ELSE 
        xml_output = .false.
      ENDIF 

    end if ! End of me==0 section

    ! Clean memory
     IF( test_img .and. me==0 ) THEN 
       call destroy_results_out(results_out_all)
       ABI_DATATYPE_DEALLOCATE(results_out_all)
     ENDIF 

    ! This synchronization is not strictly needed, but without it,
    ! there are problems with Tv1#93 in parallel, PGI compiler, on Intel/PC
    CALL abi_io_redirect(new_io_comm=xmpi_world)

    CALL timab(44,2,tsec)

    WRITE(*,*) 'End of step_12'

  END SUBROUTINE 


  SUBROUTINE test_dataset()
    
    WRITE(*,*)
    WRITE(*,*) 'size(dtsets) = ', size(dtsets)
    WRITE(*,*)
    ! dtsets is indexed by 0:ndtset_alloc


    !
    ! print out some fields that are defined in the input file
    !
    WRITE(*,*)
    WRITE(*,*) 'Dataset'
    WRITE(*,*)
    WRITE(*,*) 'natom = ', dtsets(1)%natom
    WRITE(*,*) 'ntypat = ', dtsets(1)%ntypat
    WRITE(*,*) 'znucl = ', dtsets(1)%znucl(:)
    WRITE(*,*) 'acell = ', results_out(1)%acell(:,1)
    WRITE(*,*) 'toldfe = ', dtsets(1)%toldfe

    WRITE(*,*)
    WRITE(*,*) 'Pass here ...'
  END SUBROUTINE 



END PROGRAM 

