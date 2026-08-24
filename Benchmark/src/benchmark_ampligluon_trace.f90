program benchmark_ampligluon_trace
  use ampligluon_trace, only: ampligluon_trace_t
  use benchmark_events, only: benchmark_dp, benchmark_event_t, benchmark_fail
  implicit none

  type(ampligluon_trace_t) :: amplitude
  type(benchmark_event_t), allocatable :: events(:)
  integer, allocatable :: cell_repetitions(:, :)
  logical, allocatable :: timed_cells(:, :)
  real(benchmark_dp), allocatable :: cell_seconds(:, :)
  real(benchmark_dp), allocatable :: matrix_values(:, :)
  character(len=4096) :: argument
  character(len=64) :: backend_name, colour_mode, mhv_mode
  real(benchmark_dp) :: cell_target_seconds, checksum, sample_pass_seconds
  real(benchmark_dp) :: helicity_sum_block_seconds, helicity_sum_seconds
  real(benchmark_dp) :: initialization_seconds, matrix_element
  real(benchmark_dp) :: start_time, stop_time, target_seconds
  integer :: argument_count, batches, batch, configuration, event_index
  integer :: evaluations_per_sweep, first_event_argument, ios, last_event_argument
  integer :: number_of_events, repetition_quantum
  integer :: helicity_sum_repetitions
  logical :: initialization_only, use_colour_fft, use_repetition_quantum
  logical :: sum_helicities, use_mhv_optimization

  argument_count = command_argument_count()
  if (argument_count < 5) call benchmark_fail( &
       'usage: benchmark_ampligluon_trace MHV_MODE COLOUR_MODE TARGET_SECONDS '// &
       'BATCHES EVENT... [--sum-helicities] [--repetition-quantum=N] '// &
       '[--initialization-only]')
  call get_command_argument(1, mhv_mode)
  select case (trim(mhv_mode))
  case ('optimized-mhv')
    use_mhv_optimization = .true.
  case ('default-bg')
    use_mhv_optimization = .false.
  case default
    call benchmark_fail("MHV_MODE must be 'optimized-mhv' or 'default-bg'")
  end select
  call get_command_argument(2, argument)
  colour_mode = trim(argument)
  select case (trim(colour_mode))
  case ('fft')
    use_colour_fft = .true.
    if (use_mhv_optimization) then
      backend_name = 'AmpliGluonTraceOptimizedMHV'
    else
      backend_name = 'AmpliGluonTraceDefaultBG'
    end if
  case ('direct')
    use_colour_fft = .false.
    if (use_mhv_optimization) then
      backend_name = 'AmpliGluonTraceOptimizedMHVDirectColour'
    else
      backend_name = 'AmpliGluonTraceDefaultBGDirectColour'
    end if
  case default
    call benchmark_fail("COLOUR_MODE must be 'fft' or 'direct'")
  end select
  call get_command_argument(3, argument)
  read(argument, *, iostat=ios) target_seconds
  if (ios /= 0 .or. target_seconds <= 0.0_benchmark_dp) &
       call benchmark_fail('TARGET_SECONDS must be positive')
  call get_command_argument(4, argument)
  read(argument, *, iostat=ios) batches
  if (ios /= 0 .or. batches < 1) call benchmark_fail('BATCHES must be positive')

  initialization_only = .false.
  sum_helicities = .false.
  use_repetition_quantum = .false.
  repetition_quantum = 1
  last_event_argument = argument_count
  call get_command_argument(argument_count, argument)
  if (trim(argument) == '--initialization-only') then
    initialization_only = .true.
    last_event_argument = argument_count-1
  end if
  call get_command_argument(last_event_argument, argument)
  if (index(trim(argument), '--repetition-quantum=') == 1) then
    read(argument(22:), *, iostat=ios) repetition_quantum
    if (ios /= 0 .or. repetition_quantum < 1) &
         call benchmark_fail('--repetition-quantum must be positive')
    use_repetition_quantum = .true.
    last_event_argument = last_event_argument-1
  end if
  call get_command_argument(last_event_argument, argument)
  if (trim(argument) == '--sum-helicities') then
    sum_helicities = .true.
    last_event_argument = last_event_argument-1
  end if
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
  write(*, '(a,1x,i0)') 'DIMENSION', amplitude%number_of_colour_orders()
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
           use_mhv_optimization=use_mhv_optimization)
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
  if (sum_helicities) then
    do configuration = 1, events(1)%number_of_helicities
      timed_cells(configuration, 1) = .not. analytic_zero( &
           events(1)%helicities(:, configuration))
    end do
  else
    do event_index = 1, number_of_events
      do configuration = 1, events(event_index)%number_of_helicities
        timed_cells(configuration, event_index) = .not. analytic_zero( &
             events(event_index)%helicities(:, configuration))
      end do
    end do
  end if
  if (.not. sum_helicities .and. .not. use_mhv_optimization) then
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
  if (sum_helicities) then
    if (use_repetition_quantum) then
      call calibrate_helicity_sum(repetition_quantum, &
           helicity_sum_repetitions, helicity_sum_seconds, &
           helicity_sum_block_seconds, checksum)
      call write_helicity_sum_calibration(helicity_sum_block_seconds)
      do batch = 1, batches
        call evaluate_helicity_sum(helicity_sum_repetitions, &
             helicity_sum_block_seconds, checksum)
        helicity_sum_seconds = helicity_sum_block_seconds/ &
             real(helicity_sum_repetitions, benchmark_dp)
        call write_helicity_sum_batch(batch, helicity_sum_repetitions, &
             helicity_sum_seconds)
      end do
    else
      call calibrate_helicity_sum(1, helicity_sum_repetitions, &
           helicity_sum_seconds, helicity_sum_block_seconds, checksum)
      call write_helicity_sum_calibration(helicity_sum_block_seconds)
      call write_helicity_sum_batch(1, helicity_sum_repetitions, &
           helicity_sum_seconds)
      do batch = 2, batches
        call evaluate_helicity_sum(helicity_sum_repetitions, &
             helicity_sum_block_seconds, checksum)
        helicity_sum_seconds = helicity_sum_block_seconds/ &
             real(helicity_sum_repetitions, benchmark_dp)
        call write_helicity_sum_batch(batch, helicity_sum_repetitions, &
             helicity_sum_seconds)
      end do
    end if
  else if (use_repetition_quantum) then
    call calibrate_batch(cell_target_seconds, repetition_quantum, &
         cell_repetitions, cell_seconds, checksum)
    call write_calibration(cell_repetitions, cell_seconds)
    do batch = 1, batches
      call evaluate_batch(cell_repetitions, cell_seconds, checksum)
      call write_batch(batch, cell_repetitions, cell_seconds)
    end do
  else
    call calibrate_batch(cell_target_seconds, 1, cell_repetitions, &
         cell_seconds, checksum)
    call write_batch(1, cell_repetitions, cell_seconds)
    do batch = 2, batches
      call evaluate_batch(cell_repetitions, cell_seconds, checksum)
      call write_batch(batch, cell_repetitions, cell_seconds)
    end do
  end if
  write(*, '(a,1x,es24.16)') 'CHECKSUM', checksum

contains

  logical function analytic_zero(helicities) result(is_zero)
    integer, intent(in) :: helicities(:)

    integer :: outgoing_positive

    outgoing_positive = count(helicities(1:2) < 0)+count(helicities(3:) > 0)
    is_zero = outgoing_positive < 2 .or. &
         outgoing_positive > size(helicities)-2
  end function analytic_zero

  subroutine evaluate_helicity_sum(number_of_repetitions, evaluation_seconds, &
                                   result_checksum)
    integer, intent(in) :: number_of_repetitions
    real(benchmark_dp), intent(out) :: evaluation_seconds
    real(benchmark_dp), intent(out) :: result_checksum

    real(benchmark_dp) :: after, before, summed_matrix_element
    integer :: configuration_index, repetition_index

    result_checksum = 0.0_benchmark_dp
    call cpu_time(before)
    do repetition_index = 1, number_of_repetitions
      summed_matrix_element = 0.0_benchmark_dp
      do configuration_index = 1, events(1)%number_of_helicities
        if (.not. timed_cells(configuration_index, 1)) cycle
        call amplitude%evaluate(events(1)%momenta, &
             events(1)%helicities(:, configuration_index), matrix_element, &
             strong_coupling=events(1)%strong_coupling, &
             use_mhv_optimization=use_mhv_optimization)
        summed_matrix_element = summed_matrix_element+matrix_element
      end do
      result_checksum = result_checksum+summed_matrix_element
    end do
    call cpu_time(after)
    evaluation_seconds = after-before
  end subroutine evaluate_helicity_sum

  subroutine calibrate_helicity_sum(minimum_repetitions, &
                                    number_of_repetitions, seconds_per_sum, &
                                    block_seconds, result_checksum)
    integer, intent(in) :: minimum_repetitions
    integer, intent(out) :: number_of_repetitions
    real(benchmark_dp), intent(out) :: seconds_per_sum, block_seconds
    real(benchmark_dp), intent(out) :: result_checksum

    number_of_repetitions = minimum_repetitions
    do
      call evaluate_helicity_sum(number_of_repetitions, block_seconds, &
           result_checksum)
      if (block_seconds >= target_seconds) exit
      if (number_of_repetitions > huge(number_of_repetitions)- &
           number_of_repetitions) &
           call benchmark_fail('helicity-sum repetition count overflow')
      number_of_repetitions = 2*number_of_repetitions
    end do
    seconds_per_sum = block_seconds/ &
         real(number_of_repetitions, benchmark_dp)
  end subroutine calibrate_helicity_sum

  subroutine write_helicity_sum_calibration(block_seconds)
    real(benchmark_dp), intent(in) :: block_seconds

    write(*, '(a,2(1x,i0),1x,es24.16)') &
         'CALIBRATION_CELL_TOTAL_SECONDS', 1, 1, block_seconds
  end subroutine write_helicity_sum_calibration

  subroutine write_helicity_sum_batch(batch_number, number_of_repetitions, &
                                      seconds_per_sum)
    integer, intent(in) :: batch_number, number_of_repetitions
    real(benchmark_dp), intent(in) :: seconds_per_sum

    write(*, '(a,1x,i0,1x,es24.16)') 'EVALUATION_SWEEP_SECONDS', &
         batch_number, seconds_per_sum
    write(*, '(a,3(1x,i0),1x,es24.16)') 'EVALUATION_CELL_SECONDS', &
         batch_number, 1, 1, seconds_per_sum
    if (use_repetition_quantum) &
         write(*, '(a,4(1x,i0))') 'EVALUATION_CELL_REPETITIONS', &
         batch_number, 1, 1, number_of_repetitions
  end subroutine write_helicity_sum_batch

  subroutine calibrate_batch(seconds_per_cell, minimum_repetitions, &
                             number_of_repetitions, evaluation_seconds, &
                             result_checksum)
    real(benchmark_dp), intent(in) :: seconds_per_cell
    integer, intent(in) :: minimum_repetitions
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
        number_of_repetitions(configuration, event_index) = minimum_repetitions
        do
          cell_checksum = 0.0_benchmark_dp
          call cpu_time(before)
          do repetition_index = 1, &
               number_of_repetitions(configuration, event_index)
            call amplitude%evaluate(events(event_index)%momenta, &
                 events(event_index)%helicities(:, configuration), &
                 matrix_element, &
                 strong_coupling=events(event_index)%strong_coupling, &
                 use_mhv_optimization=use_mhv_optimization)
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
               use_mhv_optimization=use_mhv_optimization)
          result_checksum = result_checksum+matrix_element
        end do
        call cpu_time(after)
        evaluation_seconds(configuration, event_index) = (after-before)/ &
             real(number_of_repetitions(configuration, event_index), benchmark_dp)
      end do
    end do
  end subroutine evaluate_batch

  subroutine write_calibration(number_of_repetitions, evaluation_seconds)
    integer, intent(in) :: number_of_repetitions(:, :)
    real(benchmark_dp), intent(in) :: evaluation_seconds(:, :)

    do event_index = 1, number_of_events
      do configuration = 1, events(event_index)%number_of_helicities
        if (.not. timed_cells(configuration, event_index)) cycle
        write(*, '(a,2(1x,i0),1x,es24.16)') &
             'CALIBRATION_CELL_TOTAL_SECONDS', event_index, configuration, &
             evaluation_seconds(configuration, event_index)* &
             real(number_of_repetitions(configuration, event_index), benchmark_dp)
      end do
    end do
  end subroutine write_calibration

  subroutine write_batch(batch_number, number_of_repetitions, evaluation_seconds)
    integer, intent(in) :: batch_number
    integer, intent(in) :: number_of_repetitions(:, :)
    real(benchmark_dp), intent(in) :: evaluation_seconds(:, :)

    write(*, '(a,1x,i0,1x,es24.16)') 'EVALUATION_SWEEP_SECONDS', &
         batch_number, sum(evaluation_seconds)
    do event_index = 1, number_of_events
      do configuration = 1, events(event_index)%number_of_helicities
        if (.not. timed_cells(configuration, event_index)) cycle
        write(*, '(a,3(1x,i0),1x,es24.16)') 'EVALUATION_CELL_SECONDS', &
             batch_number, event_index, configuration, &
             evaluation_seconds(configuration, event_index)
        if (use_repetition_quantum) &
             write(*, '(a,4(1x,i0))') 'EVALUATION_CELL_REPETITIONS', &
             batch_number, event_index, configuration, &
             number_of_repetitions(configuration, event_index)
      end do
    end do
  end subroutine write_batch

end program benchmark_ampligluon_trace
