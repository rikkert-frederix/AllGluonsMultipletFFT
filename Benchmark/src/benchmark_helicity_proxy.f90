program benchmark_helicity_proxy
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use benchmark_events, only: benchmark_dp, benchmark_event_t, benchmark_fail
  use gluon_kinematics, only: external_massless_vector
  use trace_order_recursion, only: evaluate_canonical_trace_order
  implicit none

  type(benchmark_event_t) :: event
  character(len=4096) :: event_file
  integer :: event_index, number_of_events

  number_of_events = command_argument_count()
  if (number_of_events < 1) call benchmark_fail( &
       'usage: benchmark_helicity_proxy EVENT...')

  write(*, '(a,1x,a)') 'LC_PROXY', 'CanonicalTraceOrderBG'
  write(*, '(a,1x,i0)') 'LC_EVENTS', number_of_events
  do event_index = 1, number_of_events
    call get_command_argument(event_index, event_file)
    call event%load(trim(event_file))
    call write_event_weights(event, event_index)
  end do

contains

  subroutine write_event_weights(input_event, index)
    type(benchmark_event_t), intent(in) :: input_event
    integer, intent(in) :: index

    complex(benchmark_dp) :: partial_amplitude
    complex(benchmark_dp), allocatable :: wavefunctions(:, :)
    real(benchmark_dp), allocatable :: outgoing_momenta(:, :), weights(:)
    real(benchmark_dp) :: helicity_sum, mhv_fraction, mhv_sum, weight
    integer :: configuration, external_gluon, outgoing_positive, total_gluons
    logical :: is_mhv, is_zero

    total_gluons = input_event%final_gluons+2
    if (.not. ieee_is_finite(input_event%strong_coupling)) &
         call benchmark_fail('strong coupling is not finite')
    if (.not. all(ieee_is_finite(input_event%momenta))) &
         call benchmark_fail('event momenta are not finite')
    if (any(abs(input_event%helicities) /= 1)) &
         call benchmark_fail('all gluon helicities must be -1 or +1')

    allocate(outgoing_momenta(0:3, total_gluons))
    allocate(wavefunctions(4, total_gluons))
    allocate(weights(input_event%number_of_helicities))
    outgoing_momenta(:, 1:2) = -input_event%momenta(:, 1:2)
    outgoing_momenta(:, 3:total_gluons) = &
         input_event%momenta(:, 3:total_gluons)

    weights = 0.0_benchmark_dp
    helicity_sum = 0.0_benchmark_dp
    mhv_sum = 0.0_benchmark_dp
    write(*, '(a,3(1x,i0))') 'LC_EVENT', index, total_gluons, &
         input_event%number_of_helicities
    do configuration = 1, input_event%number_of_helicities
      outgoing_positive = count( &
           input_event%helicities(1:2, configuration) < 0)+count( &
           input_event%helicities(3:total_gluons, configuration) > 0)
      is_zero = outgoing_positive < 2 .or. &
           outgoing_positive > total_gluons-2 .or. &
           abs(input_event%strong_coupling) <= 0.0_benchmark_dp
      is_mhv = outgoing_positive == 2 .or. &
           outgoing_positive == total_gluons-2

      if (.not. is_zero) then
        do external_gluon = 1, total_gluons
          call external_massless_vector( &
               outgoing_momenta(:, external_gluon), &
               input_event%helicities(external_gluon, configuration), &
               wavefunctions(:, external_gluon))
        end do
        call evaluate_canonical_trace_order( &
             outgoing_momenta(:, 1:total_gluons-1), &
             wavefunctions(:, 1:total_gluons-1), &
             wavefunctions(:, total_gluons), input_event%strong_coupling, &
             partial_amplitude)
        if (.not. ieee_is_finite(real(partial_amplitude, benchmark_dp)) .or. &
             .not. ieee_is_finite(aimag(partial_amplitude))) &
             call benchmark_fail('leading-colour proxy amplitude is not finite')
        weight = real(partial_amplitude, benchmark_dp)**2+ &
             aimag(partial_amplitude)**2
        if (.not. ieee_is_finite(weight) .or. weight < 0.0_benchmark_dp) &
             call benchmark_fail('leading-colour proxy weight is not finite')
        weights(configuration) = weight
        helicity_sum = helicity_sum+weight
        if (is_mhv) mhv_sum = mhv_sum+weight
        if (.not. ieee_is_finite(helicity_sum) .or. &
             .not. ieee_is_finite(mhv_sum)) &
             call benchmark_fail('leading-colour proxy sum is not finite')
      end if
      write(*, '(a,2(1x,i0),1x,es24.16)') 'LC_WEIGHT', index, &
           configuration, weights(configuration)
    end do

    mhv_fraction = 0.0_benchmark_dp
    if (helicity_sum > 0.0_benchmark_dp) &
         mhv_fraction = min(1.0_benchmark_dp, max(0.0_benchmark_dp, &
         mhv_sum/helicity_sum))
    write(*, '(a,1x,i0,1x,es24.16)') 'LC_HELICITY_SUM', index, helicity_sum
    write(*, '(a,1x,i0,1x,es24.16)') 'LC_MHV_FRACTION', index, mhv_fraction
  end subroutine write_event_weights

end program benchmark_helicity_proxy
