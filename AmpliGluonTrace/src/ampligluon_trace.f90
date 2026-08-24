module ampligluon_trace
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ampligluon_common, only: dp, fail, i64
  use gluon_kinematics, only: external_massless_vector
  use trace_colour_matrix, only: trace_colour_matrix_t
  use trace_current_dag, only: trace_current_dag_t
  use trace_mhv, only: trace_mhv_t
  use trace_order_recursion, only: trace_order_recursion_t
  implicit none
  private

  real(dp), parameter :: initial_colour_states = 64.0_dp
  real(dp), parameter :: momentum_tolerance = 2.0e-9_dp
  integer, parameter :: maximum_shared_current_degree = 10

  type, public :: ampligluon_trace_t
    private
    logical :: initialized = .false.
    logical :: use_colour_fft_value = .true.
    integer :: total_gluons = 0
    integer :: current_gluons = 0
    integer :: number_of_orders_value = 0
    type(trace_colour_matrix_t), allocatable :: colour
    type(trace_mhv_t), allocatable :: mhv
    type(trace_current_dag_t), allocatable :: current_dag
    type(trace_order_recursion_t), allocatable :: order_recursion
    real(dp), allocatable :: outgoing_momenta(:, :)
    complex(dp), allocatable :: external_wavefunctions(:, :)
    complex(dp), allocatable :: amplitudes(:)
  contains
    procedure, public :: evaluate => evaluate_ampligluon_trace
    procedure, public :: initialize => initialize_ampligluon_trace
    procedure, public :: number_of_colour_orders
    procedure, public :: number_of_final_gluons
  end type ampligluon_trace_t

contains

  subroutine initialize_ampligluon_trace(this, final_gluons, use_colour_fft)
    class(ampligluon_trace_t), intent(inout) :: this
    integer, intent(in) :: final_gluons
    logical, intent(in), optional :: use_colour_fft

    if (final_gluons < 2) &
         call fail('gg scattering needs at least two final gluons')
    this%initialized = .false.
    this%use_colour_fft_value = .true.
    if (present(use_colour_fft)) &
         this%use_colour_fft_value = use_colour_fft
    this%total_gluons = final_gluons+2
    this%current_gluons = this%total_gluons-1
    this%number_of_orders_value = &
         trace_order_count(this%current_gluons)

    if (allocated(this%colour)) deallocate(this%colour)
    if (allocated(this%mhv)) deallocate(this%mhv)
    if (allocated(this%current_dag)) deallocate(this%current_dag)
    if (allocated(this%order_recursion)) deallocate(this%order_recursion)

    if (allocated(this%outgoing_momenta)) deallocate(this%outgoing_momenta)
    if (allocated(this%external_wavefunctions)) &
         deallocate(this%external_wavefunctions)
    if (allocated(this%amplitudes)) deallocate(this%amplitudes)

    allocate(this%outgoing_momenta(0:3, this%total_gluons))
    allocate(this%external_wavefunctions(4, this%total_gluons))
    this%initialized = .true.
  end subroutine initialize_ampligluon_trace

  subroutine evaluate_ampligluon_trace(this, momenta, helicities, &
       matrix_element_squared, strong_coupling, average_initial_colours, &
       ordered_amplitudes, use_mhv_optimization)
    class(ampligluon_trace_t), intent(inout) :: this
    real(dp), intent(in) :: momenta(0:, :)
    integer, intent(in) :: helicities(:)
    real(dp), intent(out) :: matrix_element_squared
    real(dp), intent(in), optional :: strong_coupling
    logical, intent(in), optional :: average_initial_colours
    complex(dp), allocatable, intent(inout), optional :: ordered_amplitudes(:)
    logical, intent(in), optional :: use_mhv_optimization

    real(dp) :: amplitude_scale, coupling, squared_scale
    integer :: gluon, outgoing_positive
    logical :: average_colours, enable_mhv_optimization

    if (.not. this%initialized) &
         call fail('amplitude object is not initialized')
    call validate_event(this, momenta, helicities)

    coupling = 1.0_dp
    if (present(strong_coupling)) coupling = strong_coupling
    if (.not. ieee_is_finite(coupling)) &
         call fail('strong coupling is not finite')
    average_colours = .false.
    if (present(average_initial_colours)) &
         average_colours = average_initial_colours
    enable_mhv_optimization = .true.
    if (present(use_mhv_optimization)) &
         enable_mhv_optimization = use_mhv_optimization

    ! Pure-gluon tree amplitudes vanish when the corresponding all-outgoing
    ! helicities contain fewer than two states of either sign.  Crossing flips
    ! the helicities of the two incoming gluons.  Return the exact zero before
    ! constructing wavefunctions, currents, or the colour transform, while
    ! preserving validation and the optional ordered-amplitude interface.
    outgoing_positive = count(helicities(1:2) < 0)+ &
         count(helicities(3:this%total_gluons) > 0)
    if (abs(coupling) <= 0.0_dp .or. outgoing_positive < 2 .or. &
         outgoing_positive > this%total_gluons-2) then
      matrix_element_squared = 0.0_dp
      if (present(ordered_amplitudes)) then
        if (allocated(ordered_amplitudes)) then
          if (size(ordered_amplitudes) /= this%number_of_orders_value) &
               deallocate(ordered_amplitudes)
        end if
        if (.not. allocated(ordered_amplitudes)) &
             allocate(ordered_amplitudes(this%number_of_orders_value))
        ordered_amplitudes = cmplx(0.0_dp, 0.0_dp, dp)
      end if
      return
    end if
    amplitude_scale = coupling**(this%total_gluons-2)
    squared_scale = amplitude_scale*amplitude_scale
    if (.not. ieee_is_finite(amplitude_scale) .or. &
         .not. ieee_is_finite(squared_scale)) &
         call fail('strong coupling power is not finite')
    call ensure_amplitude_storage(this)

    this%outgoing_momenta(:, 1:2) = -momenta(:, 1:2)
    this%outgoing_momenta(:, 3:this%total_gluons) = &
         momenta(:, 3:this%total_gluons)
    do gluon = 1, this%total_gluons
      call external_massless_vector(this%outgoing_momenta(:, gluon), &
           helicities(gluon), this%external_wavefunctions(:, gluon))
    end do

    if (enable_mhv_optimization .and. (outgoing_positive == 2 .or. &
         outgoing_positive == this%total_gluons-2)) then
      call ensure_mhv_evaluator(this)
      call this%mhv%evaluate(this%outgoing_momenta, &
           this%external_wavefunctions, helicities, 1.0_dp, this%amplitudes)
    else
      call ensure_kinematic_recursion(this)
      if (allocated(this%current_dag)) then
        call this%current_dag%evaluate( &
             this%outgoing_momenta(:, 1:this%current_gluons), &
             this%external_wavefunctions(:, 1:this%current_gluons), &
             this%external_wavefunctions(:, this%total_gluons), 1.0_dp, &
             this%amplitudes)
      else if (allocated(this%order_recursion)) then
        call this%order_recursion%evaluate( &
             this%outgoing_momenta(:, 1:this%current_gluons), &
             this%external_wavefunctions(:, 1:this%current_gluons), &
             this%external_wavefunctions(:, this%total_gluons), 1.0_dp, &
             this%amplitudes)
      else
        call fail('trace kinematic recursion is not initialized')
      end if
    end if
    call ensure_colour_contraction(this)
    call this%colour%contract(this%amplitudes, matrix_element_squared)
    matrix_element_squared = squared_scale*matrix_element_squared
    if (.not. ieee_is_finite(matrix_element_squared)) &
         call fail('trace matrix element is not finite after coupling scaling')
    if (average_colours) &
         matrix_element_squared = matrix_element_squared/initial_colour_states

    if (present(ordered_amplitudes)) then
      if (allocated(ordered_amplitudes)) then
        if (size(ordered_amplitudes) /= size(this%amplitudes)) &
             deallocate(ordered_amplitudes)
      end if
      if (.not. allocated(ordered_amplitudes)) &
           allocate(ordered_amplitudes(size(this%amplitudes)))
      ordered_amplitudes = amplitude_scale*this%amplitudes
    end if
  end subroutine evaluate_ampligluon_trace

  subroutine ensure_amplitude_storage(this)
    class(ampligluon_trace_t), intent(inout) :: this

    if (.not. allocated(this%amplitudes)) &
         allocate(this%amplitudes(this%number_of_orders_value))
  end subroutine ensure_amplitude_storage

  subroutine ensure_colour_contraction(this)
    class(ampligluon_trace_t), intent(inout) :: this

    if (allocated(this%colour)) return
    allocate(this%colour)
    call this%colour%initialize(this%total_gluons, &
         use_colour_fft=this%use_colour_fft_value)
    if (this%colour%number_of_orders() /= this%number_of_orders_value) &
         call fail('trace colour order count changed during initialization')
  end subroutine ensure_colour_contraction

  subroutine ensure_mhv_evaluator(this)
    class(ampligluon_trace_t), intent(inout) :: this

    if (allocated(this%mhv)) return
    allocate(this%mhv)
    call this%mhv%initialize(this%total_gluons, &
         this%number_of_orders_value)
  end subroutine ensure_mhv_evaluator

  subroutine ensure_kinematic_recursion(this)
    class(ampligluon_trace_t), intent(inout) :: this

    if (allocated(this%current_dag) .or. allocated(this%order_recursion)) return
    if (this%current_gluons <= maximum_shared_current_degree) then
      allocate(this%current_dag)
      call this%current_dag%initialize(this%current_gluons)
      if (this%current_dag%number_of_orders() /= &
           this%number_of_orders_value) &
           call fail('trace current and colour order counts disagree')
    else
      allocate(this%order_recursion)
      call this%order_recursion%initialize(this%current_gluons, &
           this%number_of_orders_value)
      if (this%order_recursion%number_of_orders() /= &
           this%number_of_orders_value) &
           call fail('trace recursion and colour order counts disagree')
    end if
  end subroutine ensure_kinematic_recursion

  subroutine validate_event(this, momenta, helicities)
    class(ampligluon_trace_t), intent(in) :: this
    real(dp), intent(in) :: momenta(0:, :)
    integer, intent(in) :: helicities(:)

    real(dp) :: mass_squared, momentum_scale, residual(0:3)
    integer :: external_gluon

    if (size(momenta, 1) /= 4 .or. &
         size(momenta, 2) /= this%total_gluons) &
         call fail('momentum array has the wrong shape')
    if (size(helicities) /= this%total_gluons) &
         call fail('helicity array has the wrong length')
    if (.not. all(ieee_is_finite(momenta))) &
         call fail('momenta are not finite')
    if (any(abs(helicities) /= 1)) &
         call fail('all gluon helicities must be -1 or +1')

    momentum_scale = max(1.0_dp, maxval(abs(momenta)))
    do external_gluon = 1, this%total_gluons
      if (momenta(0, external_gluon) <= 0.0_dp) &
           call fail('physical external energies must be positive')
      mass_squared = momenta(0, external_gluon)**2- &
           sum(momenta(1:3, external_gluon)**2)
      if (abs(mass_squared) > momentum_tolerance*momentum_scale**2) &
           call fail('external momentum is not massless')
    end do
    residual = momenta(:, 1)+momenta(:, 2)- &
         sum(momenta(:, 3:this%total_gluons), dim=2)
    if (maxval(abs(residual)) > momentum_tolerance*momentum_scale) &
         call fail('external momenta do not conserve four-momentum')
  end subroutine validate_event

  integer function number_of_colour_orders(this) result(number)
    class(ampligluon_trace_t), intent(in) :: this

    if (.not. this%initialized) &
         call fail('amplitude object is not initialized')
    number = this%number_of_orders_value
  end function number_of_colour_orders

  integer function number_of_final_gluons(this) result(number)
    class(ampligluon_trace_t), intent(in) :: this

    if (.not. this%initialized) &
         call fail('amplitude object is not initialized')
    number = this%total_gluons-2
  end function number_of_final_gluons

  integer function trace_order_count(number) result(value)
    integer, intent(in) :: number

    integer(i64) :: value_64
    integer :: factor

    if (number < 0) call fail('negative trace-order factorial')
    value_64 = 1_i64
    do factor = 2, number
      if (value_64 > huge(value_64)/int(factor, i64)) &
           call fail('the number of colour orders overflows a 64-bit integer')
      value_64 = value_64*int(factor, i64)
    end do
    if (value_64 > int(huge(value), i64)) &
         call fail('the number of colour orders exceeds array-index capacity')
    value = int(value_64)
  end function trace_order_count

end module ampligluon_trace
