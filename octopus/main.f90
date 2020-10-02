! Various macros

!-----------------------------------------------------------------------
program octopus_debug
!-----------------------------------------------------------------------

  use utils_oct_m
  use command_line_oct_m
  use block_t_oct_m
  use global_oct_m
  use messages_oct_m
  use varinfo_oct_m
  use parser_oct_m
  use run_oct_m
  use io_oct_m
  use calc_mode_par_oct_m
  use profiling_oct_m

  implicit none

  ! Local variables
  character(len=256) :: config_str
  integer :: inp_calc_mode, ierr
  type(block_t) :: blk

  call getopt_init(ierr)
  config_str = trim(get_config_opts()) // trim(get_optional_libraries())
  if(ierr  ==  0) call getopt_octopus(trim(config_str))
  call getopt_end()

  write(*,*) 'config_str = ', trim(config_str)

  call global_init()
  call messages_init()

  call parse_variable('ReportMemory', .false., conf%report_memory)

  if(parse_block('CalculationMode', blk) == 0) then
    call messages_write('The datasets mode has been deprecated,', new_line = .true.)
    call messages_write('please use several Octopus runs.')
    call messages_fatal()
  end if

  call parse_variable('CalculationMode', CM_GS, inp_calc_mode)
  if(.not.varinfo_valid_option('CalculationMode', inp_calc_mode)) call messages_input_error('CalculationMode')

  ! Now we can initialize the I/O
  call io_init()

  call calc_mode_par_init()

  ! now we declare octopus as running
  call messages_switch_status('running')

  call profiling_init()

  call print_header()

  write(*,*) 'inp_calc_mode = ', inp_calc_mode

  !call run(inp_calc_mode)

  ! run finished successfully
  !call messages_switch_status('finished')
  !call io_end()


end program