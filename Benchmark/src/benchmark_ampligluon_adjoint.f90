program benchmark_ampligluon_adjoint
  use ampligluon_adjoint, only: ampligluon_adjoint_t
  use benchmark_events, only: benchmark_dp, benchmark_event_t, benchmark_fail
  implicit none

  type(ampligluon_adjoint_t) :: amplitude
  type(benchmark_event_t), allocatable :: events(:)
  integer, allocatable :: cell_repetitions(:, :)
  logical, allocatable :: timed_cells(:, :)
  real(benchmark_dp), allocatable :: cell_seconds(:, :)
  real(benchmark_dp), allocatable :: matrix_values(:, :)
  character(len=4096) :: argument
  real(benchmark_dp) :: cell_target_seconds, checksum, sample_pass_seconds
  real(benchmark_dp) :: initialization_seconds, matrix_element
  real(benchmark_dp) :: start_time, stop_time, target_seconds
  integer :: argument_count, batches, batch, configuration, event_index
  integer :: evaluations_per_sweep, first_event_argument, ios, last_event_argument
  integer :: number_of_events
  logical :: initialization_only, use_analytic_mhv, use_colour_fft
  character(len=64) :: backend_name, colour_mode, mhv_mode

  argument_count = command_argument_count()
  if (argument_count < 5) call benchmark_fail( &
       'usage: benchmark_ampligluon_adjoint MHV_MODE COLOUR_MODE '// &
       'TARGET_SECONDS BATCHES EVENT... [--initialization-only]')
  call get_command_argument(1, mhv_mode)
  select case (trim(mhv_mode))
  case ('optimized-mhv')
    use_analytic_mhv = .true.
  case ('default-bg')
    use_analytic_mhv = .false.
  case default
    call benchmark_fail('MHV_MODE must be optimized-mhv or default-bg')
  end select
  call get_command_argument(2, argument)
  colour_mode = trim(argument)
  select case (trim(colour_mode))
  case ('fft')
    use_colour_fft = .true.
    if (use_analytic_mhv) then
      backend_name = 'AmpliGluonAdjointOptimizedMHV'
    else
      backend_name = 'AmpliGluonAdjointDefaultBG'
    end if
  case ('direct')
    use_colour_fft = .false.
    if (use_analytic_mhv) then
      backend_name = 'AmpliGluonAdjointOptimizedMHVDirectColour'
    else
      backend_name = 'AmpliGluonAdjointDefaultBGDirectColour'
    end if
  case default
    call benchmark_fail('COLOUR_MODE must be fft or direct')
  end select
  call get_command_argument(3, argument)
  read(argument, *, iostat=ios) target_seconds
  if (ios /= 0 .or. target_seconds <= 0.0_benchmark_dp) &
       call benchmark_fail('TARGET_SECONDS must be positive')
  call get_command_argument(4, argument)
  read(argument, *, iostat=ios) batches
  if (ios /= 0 .or. batches < 1) call benchmark_fail('BATCHES must be positive')

  initialization_only = .false.
  last_event_argument = argument_count
  do while (last_event_argument >= 5)
    call get_command_argument(last_event_argument, argument)
    select case (trim(argument))
    case ('--initialization-only')
      if (initialization_only) &
           call benchmark_fail('duplicate --initialization-only option')
      initialization_only = .true.
    case default
      exit
    end select
    last_event_argument = last_event_argument-1
  end do
  first_event_argument = 5
  number_of_events = last_event_argument-first_event_argument+1
  if (number_of_events < 1) call benchmark_fail('at least one EVENT is required')
  allocate(events(number_of_events))
  do event_index = 1, number_of_events
    call get_command_argument(first_event_argument+event_index-1, argument)
    call events(event_index)%load(trim(argument))
    if (event_index > 1) then
      if (events(event_index)%final_gluons /= events(1)%final_gluons .or. &
           events(event_index)%number_of_helicities /= &
           events(1)%number_of_helicities) &
           call benchmark_fail('all events must have the same dimensions')
    end if
  end do

  call cpu_time(start_time)
  call amplitude%initialize(events(1)%final_gluons, &
       use_colour_fft=use_colour_fft)
  call cpu_time(stop_time)
  initialization_seconds = stop_time-start_time
  write(*, '(a,1x,a)') 'BACKEND', trim(backend_name)
  write(*, '(a,1x,i0)') 'TOTAL_GLUONS', events(1)%final_gluons+2
  write(*, '(a,1x,i0)') 'DIMENSION', amplitude%number_of_basis_amplitudes()
  write(*, '(a,1x,es24.16)') 'INITIALIZATION_SECONDS', initialization_seconds
  if (initialization_only) stop

  allocate(matrix_values(events(1)%number_of_helicities, number_of_events))
  call cpu_time(start_time)
  do event_index = 1, number_of_events
    do configuration = 1, events(event_index)%number_of_helicities
      call amplitude%evaluate(events(event_index)%momenta, &
           events(event_index)%helicities(:, configuration), &
           matrix_values(configuration, event_index), &
           strong_coupling=events(event_index)%strong_coupling, &
           use_analytic_mhv=use_analytic_mhv)
    end do
  end do
  call cpu_time(stop_time)
  sample_pass_seconds = stop_time-start_time
  write(*, '(a,1x,es24.16)') 'FIRST_SAMPLE_PASS_SECONDS', &
       sample_pass_seconds
  do event_index = 1, number_of_events
    do configuration = 1, events(event_index)%number_of_helicities
      write(*, '(a,2(1x,i0),1x,es24.16)') 'MATRIX_ELEMENT', &
           event_index, configuration, matrix_values(configuration, event_index)
    end do
  end do

  allocate(timed_cells(events(1)%number_of_helicities, number_of_events))
  timed_cells = .false.
  do event_index = 1, number_of_events
    do configuration = 1, events(event_index)%number_of_helicities
      timed_cells(configuration, event_index) = .not. analytic_zero( &
           events(event_index)%helicities(:, configuration))
    end do
  end do
  if (.not. use_analytic_mhv) then
    find_event: do event_index = 1, number_of_events
      do configuration = 1, events(event_index)%number_of_helicities
        if (timed_cells(configuration, event_index)) then
          timed_cells = .false.
          timed_cells(configuration, event_index) = .true.
          exit find_event
        end if
      end do
    end do find_event
  end if
  evaluations_per_sweep = count(timed_cells)
  if (evaluations_per_sweep == 0) &
       call benchmark_fail('sample contains no nonzero matrix element')
  write(*, '(a,1x,i0)') 'EVALUATIONS_PER_SWEEP', evaluations_per_sweep
  cell_target_seconds = target_seconds/real(evaluations_per_sweep, benchmark_dp)
  allocate(cell_repetitions(events(1)%number_of_helicities, number_of_events))
  allocate(cell_seconds(events(1)%number_of_helicities, number_of_events))
  call calibrate_batch(cell_target_seconds, cell_repetitions, cell_seconds, &
       checksum)
  call write_batch(1, cell_seconds)
  do batch = 2, batches
    call evaluate_batch(cell_repetitions, cell_seconds, checksum)
    call write_batch(batch, cell_seconds)
  end do
  write(*, '(a,1x,es24.16)') 'CHECKSUM', checksum

contains

  logical function analytic_zero(helicities) result(is_zero)
    integer, intent(in) :: helicities(:)

    integer :: outgoing_positive

    outgoing_positive = count(helicities(1:2) < 0)+count(helicities(3:) > 0)
    is_zero = outgoing_positive < 2 .or. &
         outgoing_positive > size(helicities)-2
  end function analytic_zero

  subroutine calibrate_batch(seconds_per_cell, number_of_repetitions, &
                             evaluation_seconds, result_checksum)
    real(benchmark_dp), intent(in) :: seconds_per_cell
    integer, intent(out) :: number_of_repetitions(:, :)
    real(benchmark_dp), intent(out) :: evaluation_seconds(:, :)
    real(benchmark_dp), intent(out) :: result_checksum

    real(benchmark_dp) :: after, before, block_seconds, cell_checksum
    integer :: repetition_index

    result_checksum = 0.0_benchmark_dp
    number_of_repetitions = 0
    evaluation_seconds = 0.0_benchmark_dp
    do event_index = 1, number_of_events
      do configuration = 1, events(event_index)%number_of_helicities
        if (.not. timed_cells(configuration, event_index)) cycle
        number_of_repetitions(configuration, event_index) = 1
        do
          cell_checksum = 0.0_benchmark_dp
          call cpu_time(before)
          do repetition_index = 1, &
               number_of_repetitions(configuration, event_index)
            call amplitude%evaluate(events(event_index)%momenta, &
                 events(event_index)%helicities(:, configuration), &
                 matrix_element, &
                 strong_coupling=events(event_index)%strong_coupling, &
                 use_analytic_mhv=use_analytic_mhv)
            cell_checksum = cell_checksum+matrix_element
          end do
          call cpu_time(after)
          block_seconds = after-before
          if (block_seconds >= seconds_per_cell) exit
          if (number_of_repetitions(configuration, event_index) > &
               huge(number_of_repetitions(configuration, event_index))- &
               number_of_repetitions(configuration, event_index)) &
               call benchmark_fail('evaluation repetition count overflow')
          number_of_repetitions(configuration, event_index) = &
               2*number_of_repetitions(configuration, event_index)
        end do
        evaluation_seconds(configuration, event_index) = block_seconds/ &
             real(number_of_repetitions(configuration, event_index), benchmark_dp)
        result_checksum = result_checksum+cell_checksum
      end do
    end do
  end subroutine calibrate_batch

  subroutine evaluate_batch(number_of_repetitions, evaluation_seconds, &
                            result_checksum)
    integer, intent(in) :: number_of_repetitions(:, :)
    real(benchmark_dp), intent(out) :: evaluation_seconds(:, :)
    real(benchmark_dp), intent(out) :: result_checksum

    real(benchmark_dp) :: after, before
    integer :: repetition_index

    result_checksum = 0.0_benchmark_dp
    evaluation_seconds = 0.0_benchmark_dp
    do event_index = 1, number_of_events
      do configuration = 1, events(event_index)%number_of_helicities
        if (.not. timed_cells(configuration, event_index)) cycle
        call cpu_time(before)
        do repetition_index = 1, &
             number_of_repetitions(configuration, event_index)
          call amplitude%evaluate(events(event_index)%momenta, &
               events(event_index)%helicities(:, configuration), matrix_element, &
               strong_coupling=events(event_index)%strong_coupling, &
               use_analytic_mhv=use_analytic_mhv)
          result_checksum = result_checksum+matrix_element
        end do
        call cpu_time(after)
        evaluation_seconds(configuration, event_index) = (after-before)/ &
             real(number_of_repetitions(configuration, event_index), benchmark_dp)
      end do
    end do
  end subroutine evaluate_batch

  subroutine write_batch(batch_number, evaluation_seconds)
    integer, intent(in) :: batch_number
    real(benchmark_dp), intent(in) :: evaluation_seconds(:, :)

    write(*, '(a,1x,i0,1x,es24.16)') 'EVALUATION_SWEEP_SECONDS', &
         batch_number, sum(evaluation_seconds)
    do event_index = 1, number_of_events
      do configuration = 1, events(event_index)%number_of_helicities
        if (.not. timed_cells(configuration, event_index)) cycle
        write(*, '(a,3(1x,i0),1x,es24.16)') 'EVALUATION_CELL_SECONDS', &
             batch_number, event_index, configuration, &
             evaluation_seconds(configuration, event_index)
      end do
    end do
  end subroutine write_batch

end program benchmark_ampligluon_adjoint
