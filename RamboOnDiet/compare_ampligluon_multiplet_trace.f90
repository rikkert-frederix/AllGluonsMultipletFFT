program compare_ampligluon_multiplet_trace
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: error_unit
  use ampligluon_multiplet, only: ampligluon_multiplet_t
  use ampligluon_multiplet_kinds, only: dp
  use ampligluon_trace, only: ampligluon_trace_t
  use event_input, only: event_input_t
  implicit none

  real(dp), parameter :: relative_tolerance = 1.0e-10_dp
  real(dp), parameter :: absolute_zero_tolerance = 1.0e-24_dp

  type(ampligluon_multiplet_t) :: multiplet_amplitude
  type(ampligluon_trace_t) :: trace_amplitude
  type(event_input_t) :: event
  character(len=4096) :: event_file, table_file
  complex(dp), allocatable :: basis_amplitudes(:)
  integer :: configuration
  real(dp) :: difference, multiplet_result, trace_result
  real(dp) :: relative_difference, scale

  if (command_argument_count() /= 2) then
    write(error_unit, '(a)') &
         'usage: compare_ampligluon_multiplet_trace WIGNER_TABLE EVENT_FILE'
    error stop 2
  end if
  call get_command_argument(1, table_file)
  call get_command_argument(2, event_file)

  call event%load(trim(event_file))
  call multiplet_amplitude%initialize(event%final_gluons, trim(table_file))
  call trace_amplitude%initialize(event%final_gluons)

  write(*, '(a,1x,a)') 'EVENT', trim(event_file)
  do configuration = 1, event%number_of_helicities
    call multiplet_amplitude%evaluate(event%momenta, &
         event%helicities(:, configuration), basis_amplitudes, &
         multiplet_result, strong_coupling=event%strong_coupling, &
         use_mhv_optimization=.false.)
    call trace_amplitude%evaluate(event%momenta, &
         event%helicities(:, configuration), trace_result, &
         strong_coupling=event%strong_coupling, &
         use_mhv_optimization=.false.)

    if (.not. ieee_is_finite(multiplet_result) .or. &
         .not. ieee_is_finite(trace_result)) &
         call fail('a matrix element is not finite')
    difference = abs(multiplet_result-trace_result)
    scale = max(abs(multiplet_result), abs(trace_result))
    if (scale > absolute_zero_tolerance) then
      relative_difference = difference/scale
      if (relative_difference > relative_tolerance) &
           call fail('matrix elements do not agree')
    else
      relative_difference = 0.0_dp
      if (difference > absolute_zero_tolerance) &
           call fail('nominally zero matrix elements do not agree')
    end if

    write(*, '(a,1x,i0)') 'CONFIGURATION', configuration
    write(*, '(a,*(1x,i0))') 'HELICITIES', &
         event%helicities(:, configuration)
    write(*, '(a,1x,es24.16)') 'AMPLIGLUON_MULTIPLET', multiplet_result
    write(*, '(a,1x,es24.16)') 'AMPLIGLUON_TRACE', trace_result
    write(*, '(a,1x,es12.4)') 'ABSOLUTE_DIFFERENCE', difference
    write(*, '(a,1x,es12.4)') 'RELATIVE_DIFFERENCE', relative_difference
  end do
  write(*, '(a)') 'AmpliGluonMultiplet--AmpliGluonTrace comparison: PASS'

contains

  subroutine fail(message)
    character(len=*), intent(in) :: message

    write(error_unit, '(a)') &
         'compare_ampligluon_multiplet_trace: '//trim(message)
    error stop 1
  end subroutine fail

end program compare_ampligluon_multiplet_trace
