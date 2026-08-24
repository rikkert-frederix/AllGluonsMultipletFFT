program benchmark_madgraph_fixed_helicity
  use benchmark_events, only: benchmark_dp, benchmark_event_t, benchmark_fail
  implicit none

  type(benchmark_event_t), allocatable :: events(:)
  character(len=4096) :: argument
  integer, allocatable :: colours(:), repetitions(:)
  real(benchmark_dp), allocatable :: seconds(:), values(:)
  real(benchmark_dp) :: checksum, initialization_seconds, start_time
  real(benchmark_dp) :: stop_time, target_seconds
  integer :: argument_count, batch, batches, event_index, ios
  integer :: number_of_events
  real(benchmark_dp), external :: matrix

  interface
    subroutine benchmark_set_strong_coupling(strong_coupling)
      import benchmark_dp
      real(benchmark_dp), intent(in) :: strong_coupling
    end subroutine benchmark_set_strong_coupling
  end interface

  argument_count = command_argument_count()
  if (argument_count < 3) call benchmark_fail( &
       'usage: benchmark_madgraph_fixed_helicity TARGET_SECONDS BATCHES EVENT...')
  call get_command_argument(1, argument)
  read(argument, *, iostat=ios) target_seconds
  if (ios /= 0 .or. target_seconds <= 0.0_benchmark_dp) &
       call benchmark_fail('TARGET_SECONDS must be positive')
  call get_command_argument(2, argument)
  read(argument, *, iostat=ios) batches
  if (ios /= 0 .or. batches < 1) &
       call benchmark_fail('BATCHES must be positive')
  number_of_events = argument_count-2
  allocate(events(number_of_events))
  do event_index = 1, number_of_events
    call get_command_argument(event_index+2, argument)
    call events(event_index)%load(trim(argument))
    if (events(event_index)%number_of_helicities /= 1) &
         call benchmark_fail('MadGraph timing events must contain one helicity')
    if (event_index > 1) then
      if (events(event_index)%final_gluons /= events(1)%final_gluons) &
           call benchmark_fail('all events must have the same multiplicity')
      if (abs(events(event_index)%strong_coupling- &
           events(1)%strong_coupling) > epsilon(1.0_benchmark_dp)* &
           max(1.0_benchmark_dp, abs(events(1)%strong_coupling))) &
           call benchmark_fail('all events must have the same strong coupling')
    end if
    if (any(abs(events(event_index)%helicities(:, 1)) /= 1)) &
         call benchmark_fail('MadGraph helicities must be -1 or +1')
  end do

  allocate(colours(events(1)%final_gluons+2))
  colours = 1
  call cpu_time(start_time)
  call benchmark_set_strong_coupling(events(1)%strong_coupling)
  call cpu_time(stop_time)
  initialization_seconds = stop_time-start_time

  write(*, '(a,1x,a)') 'BACKEND', 'MadGraph5_aMCatNLOFixedHelicity'
  write(*, '(a,1x,a)') 'HELICITY_EVALUATOR', 'MATRIX_DIRECT_VECTOR'
  write(*, '(a,1x,i0)') 'TOTAL_GLUONS', events(1)%final_gluons+2
  write(*, '(a,1x,es24.16)') 'INITIALIZATION_SECONDS', initialization_seconds

  allocate(values(number_of_events))
  call cpu_time(start_time)
  do event_index = 1, number_of_events
    values(event_index) = matrix(events(event_index)%momenta, &
         events(event_index)%helicities(:, 1), colours)
  end do
  call cpu_time(stop_time)
  write(*, '(a,1x,es24.16)') 'FIRST_SAMPLE_PASS_SECONDS', stop_time-start_time
  do event_index = 1, number_of_events
    write(*, '(a,2(1x,i0),1x,es24.16)') 'MATRIX_ELEMENT', &
         event_index, 1, values(event_index)
  end do

  write(*, '(a,1x,i0)') 'EVALUATIONS_PER_SWEEP', number_of_events
  allocate(repetitions(number_of_events), seconds(number_of_events))
  call calibrate_batch(target_seconds/real(number_of_events, benchmark_dp), &
       repetitions, seconds, checksum)
  call write_batch(1, seconds)
  do batch = 2, batches
    call evaluate_batch(repetitions, seconds, checksum)
    call write_batch(batch, seconds)
  end do
  write(*, '(a,1x,es24.16)') 'CHECKSUM', checksum

contains

  subroutine calibrate_batch(seconds_per_cell, number_of_repetitions, &
                             evaluation_seconds, result_checksum)
    real(benchmark_dp), intent(in) :: seconds_per_cell
    integer, intent(out) :: number_of_repetitions(:)
    real(benchmark_dp), intent(out) :: evaluation_seconds(:)
    real(benchmark_dp), intent(out) :: result_checksum

    real(benchmark_dp) :: after, before, block_seconds, cell_checksum
    integer :: repetition_index

    result_checksum = 0.0_benchmark_dp
    do event_index = 1, number_of_events
      number_of_repetitions(event_index) = 1
      do
        cell_checksum = 0.0_benchmark_dp
        call cpu_time(before)
        do repetition_index = 1, number_of_repetitions(event_index)
          cell_checksum = cell_checksum+matrix(events(event_index)%momenta, &
               events(event_index)%helicities(:, 1), colours)
        end do
        call cpu_time(after)
        block_seconds = after-before
        if (block_seconds >= seconds_per_cell) exit
        if (number_of_repetitions(event_index) > &
             huge(number_of_repetitions(event_index))- &
             number_of_repetitions(event_index)) &
             call benchmark_fail('evaluation repetition count overflow')
        number_of_repetitions(event_index) = &
             2*number_of_repetitions(event_index)
      end do
      evaluation_seconds(event_index) = block_seconds/ &
           real(number_of_repetitions(event_index), benchmark_dp)
      result_checksum = result_checksum+cell_checksum
    end do
  end subroutine calibrate_batch

  subroutine evaluate_batch(number_of_repetitions, evaluation_seconds, &
                            result_checksum)
    integer, intent(in) :: number_of_repetitions(:)
    real(benchmark_dp), intent(out) :: evaluation_seconds(:)
    real(benchmark_dp), intent(out) :: result_checksum

    real(benchmark_dp) :: after, before
    integer :: repetition_index

    result_checksum = 0.0_benchmark_dp
    do event_index = 1, number_of_events
      call cpu_time(before)
      do repetition_index = 1, number_of_repetitions(event_index)
        result_checksum = result_checksum+matrix(events(event_index)%momenta, &
             events(event_index)%helicities(:, 1), colours)
      end do
      call cpu_time(after)
      evaluation_seconds(event_index) = (after-before)/ &
           real(number_of_repetitions(event_index), benchmark_dp)
    end do
  end subroutine evaluate_batch

  subroutine write_batch(batch_number, evaluation_seconds)
    integer, intent(in) :: batch_number
    real(benchmark_dp), intent(in) :: evaluation_seconds(:)

    write(*, '(a,1x,i0,1x,es24.16)') 'EVALUATION_SWEEP_SECONDS', &
         batch_number, sum(evaluation_seconds)
    do event_index = 1, number_of_events
      write(*, '(a,3(1x,i0),1x,es24.16)') 'EVALUATION_CELL_SECONDS', &
           batch_number, event_index, 1, evaluation_seconds(event_index)
    end do
  end subroutine write_batch

end program benchmark_madgraph_fixed_helicity
