#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

PROGRAM abinit
  
  USE m_argparse, ONLY : args_t
  USE m_argparse, ONLY : args_parser

  USE defs_basis, ONLY : dp, fnlen, ab_out, std_out, ch10, strlen
  USE defs_basis, ONLY : abi_io_redirect
  
  USE m_xmpi, ONLY : xmpi_world
  USE m_xmpi, ONLY : xmpi_init, xmpi_comm_rank, xmpi_end

  IMPLICIT NONE 

  INTEGER :: me
  TYPE(args_t) :: args
  CHARACTER(len=fnlen) :: filnam(5), filstat(5)
  CHARACTER(len=5000) :: message
  REAL(DP) :: tsec(2), tcpui, twalli
  CHARACTER(len=24) :: start_datetime
  INTEGER :: ndtset
  CHARACTER(len=strlen) :: string
  INTEGER :: lenstr

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


  IF( me == 0 ) CLOSE(unit=ab_out)

100 CALL xmpi_end()


CONTAINS 


  SUBROUTINE step_02()

    USE m_xpapi, ONLY : xpapi_init, xpapi_show_info
    USE m_time, ONLY : asctime

    CALL xpapi_init()
    CALL xpapi_show_info(unit=std_out,mode_paral="COLL")

    start_datetime = asctime()

    CALL timein(tcpui,twalli)

    CALL timab(1,0,tsec)

    !Start to accumulate time for the entire run. The end of accumulation is in timana.f
    CALL timab(1,1,tsec)

  END SUBROUTINE 


  SUBROUTINE step_03()
    CALL timab( 41, 3, tsec )
    ! Greetings to interactive user
    CALL iofn1( filnam, filstat, xmpi_world )
  END SUBROUTINE 


  SUBROUTINE step_04()

    USE m_build_info, ONLY : abinit_version
    USE m_build_info, ONLY : dump_config
    USE m_optim_dumper, ONLY : dump_optim
    USE m_cppopts_dumper, ONLY : dump_cpp_options
    USE m_nctk, ONLY : nctk_test_mpiio
    
    INTEGER :: ios
    CHARACTER(len=24) :: codename

    IF( me==0 ) THEN 
      OPEN( unit=ab_out, file=filnam(2),form='formatted',status='new', &
            action='write', iomsg=message, iostat=ios)
      
      IF( ios /= 0 ) THEN 
        WRITE(*,*) 'IO error:, ios = ', ios
        WRITE(*,*) 'message = ', message
        CALL xmpi_end()
      ENDIF 

      REWIND( unit=ab_out )

      codename = 'ABINIT'//repeat(' ',18)
      
      CALL herald( codename, abinit_version, ab_out )
      CALL herald( codename, abinit_version, std_out )
      
      CALL dump_config( std_out )
      CALL dump_optim( std_out )
      CALL dump_cpp_options( std_out )
      
      !
      ! Write names of files
      !
      WRITE(message, '(a,a,a,a,a,a,a,a,a,a,a,a)' ) &
          &   '- input  file    -> ',trim(filnam(1)),ch10,&
          &   '- output file    -> ',trim(filnam(2)),ch10,&
          &   '- root for input  files -> ',trim(filnam(3)),ch10,&
          &   '- root for output files -> ',trim(filnam(4)),ch10
      CALL wrtout( ab_out, message, 'COLL' )
      CALL wrtout( std_out, message, 'COLL' )
    ENDIF 

    ! Test if the netcdf library supports MPI-IO
    CALL nctk_test_mpiio()

  END SUBROUTINE 



  !5) Read the file, stringify it and return the number of datasets.
  SUBROUTINE step_05()
    CALL parsefile( filnam(1), lenstr, ndtset, string, xmpi_world )
  END SUBROUTINE 

  
  SUBROUTINE step_06_11()
    
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

    CALL memory_eval(dtsets,ab_out,mpi_enregs,ndtset,ndtset_alloc,npsp,pspheads)
  END SUBROUTINE 


END PROGRAM 

