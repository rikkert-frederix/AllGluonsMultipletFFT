module ampligluon_adjoint
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use adjoint_current_dag, only: adjoint_current_dag_t
  use adjoint_colour_matrix, only: adjoint_colour_matrix_t
  use ampligluon_common, only: dp, fail
  use gluon_kinematics, only: aux_tensor_gluon_to_gluon, &
       external_massless_vector, gluon_aux_tensor_to_gluon, &
       massless_vector_propagator, terminal_vector_contract, three_gluon, &
       two_gluon_to_aux_tensor
  implicit none
  private

  real(dp), parameter :: initial_colour_states = 64.0_dp
  real(dp), parameter :: momentum_tolerance = 2.0e-9_dp
  integer, parameter :: maximum_shared_current_degree = 10

  type, public :: ampligluon_adjoint_t
    private
    logical :: initialized = .false.
    integer :: total_gluons = 0
    integer :: current_gluons = 0
    type(adjoint_colour_matrix_t) :: colour
    type(adjoint_current_dag_t), allocatable :: current_dag
    integer, allocatable :: permutations(:, :)
    real(dp), allocatable :: outgoing_momenta(:, :)
    real(dp), allocatable :: interval_momenta(:, :, :)
    complex(dp), allocatable :: external_wavefunctions(:, :)
    complex(dp), allocatable :: gluon_current(:, :, :)
    complex(dp), allocatable :: tensor_current(:, :, :)
    complex(dp), allocatable :: amplitudes(:)
  contains
    procedure, public :: evaluate => evaluate_ampligluon_adjoint
    procedure, public :: initialize => initialize_ampligluon_adjoint
    procedure, public :: number_of_basis_amplitudes
    procedure, public :: number_of_final_gluons
  end type ampligluon_adjoint_t

contains

  subroutine initialize_ampligluon_adjoint(this, final_gluons, use_colour_fft)
    class(ampligluon_adjoint_t), intent(inout) :: this
    integer, intent(in) :: final_gluons
    logical, intent(in), optional :: use_colour_fft

    integer :: number_of_basis_elements
    logical :: enable_colour_fft

    if (final_gluons < 2) &
         call fail('gg scattering needs at least two final gluons')
    this%initialized = .false.
    this%total_gluons = final_gluons+2
    this%current_gluons = this%total_gluons-1

    enable_colour_fft = .true.
    if (present(use_colour_fft)) enable_colour_fft = use_colour_fft
    call this%colour%initialize(this%total_gluons, this%permutations, &
         use_colour_fft=enable_colour_fft)
    number_of_basis_elements = this%colour%number_of_basis_elements()

    if (allocated(this%current_dag)) deallocate(this%current_dag)

    if (allocated(this%outgoing_momenta)) deallocate(this%outgoing_momenta)
    if (allocated(this%interval_momenta)) deallocate(this%interval_momenta)
    if (allocated(this%external_wavefunctions)) &
         deallocate(this%external_wavefunctions)
    if (allocated(this%gluon_current)) deallocate(this%gluon_current)
    if (allocated(this%tensor_current)) deallocate(this%tensor_current)
    if (allocated(this%amplitudes)) deallocate(this%amplitudes)

    allocate(this%outgoing_momenta(0:3, this%total_gluons))
    allocate(this%interval_momenta(0:3, this%current_gluons, &
         this%current_gluons))
    allocate(this%external_wavefunctions(4, this%total_gluons))
    allocate(this%gluon_current(4, this%current_gluons, &
         this%current_gluons))
    allocate(this%tensor_current(6, this%current_gluons, &
         this%current_gluons))
    allocate(this%amplitudes(number_of_basis_elements))
    this%initialized = .true.
  end subroutine initialize_ampligluon_adjoint

  subroutine evaluate_ampligluon_adjoint(this, momenta, helicities, &
       matrix_element_squared, strong_coupling, average_initial_colours, &
       adjoint_amplitudes, use_analytic_mhv)
    class(ampligluon_adjoint_t), intent(inout) :: this
    real(dp), intent(in) :: momenta(0:, :)
    integer, intent(in) :: helicities(:)
    real(dp), intent(out) :: matrix_element_squared
    real(dp), intent(in), optional :: strong_coupling
    logical, intent(in), optional :: average_initial_colours
    complex(dp), allocatable, intent(inout), optional :: adjoint_amplitudes(:)
    logical, intent(in), optional :: use_analytic_mhv

    real(dp) :: coupling
    integer :: gluon, order, outgoing_positive
    logical :: analytic_mhv_enabled, average_colours, zero_coupling

    if (.not. this%initialized) &
         call fail('amplitude object is not initialized')
    call validate_event(this, momenta, helicities)

    coupling = 1.0_dp
    if (present(strong_coupling)) coupling = strong_coupling
    if (.not. ieee_is_finite(coupling)) &
         call fail('strong coupling is not finite')
    zero_coupling = .not. (coupling < 0.0_dp .or. coupling > 0.0_dp)
    average_colours = .false.
    if (present(average_initial_colours)) &
         average_colours = average_initial_colours
    analytic_mhv_enabled = .true.
    if (present(use_analytic_mhv)) analytic_mhv_enabled = use_analytic_mhv

    ! Pure-gluon tree amplitudes vanish when their all-outgoing helicities
    ! contain fewer than two states of either sign.  Crossing flips the two
    ! incoming helicities.  Avoid the high-multiplicity current DAG and colour
    ! transform entirely for these exact zero sectors and for zero coupling,
    ! while retaining the documented optional adjoint-amplitude output.
    outgoing_positive = count(helicities(1:2) < 0)+ &
         count(helicities(3:this%total_gluons) > 0)
    if (outgoing_positive < 2 .or. &
         outgoing_positive > this%total_gluons-2 .or. zero_coupling) then
      this%amplitudes = cmplx(0.0_dp, 0.0_dp, dp)
      matrix_element_squared = 0.0_dp
      if (present(adjoint_amplitudes)) &
           call copy_adjoint_amplitudes(this, adjoint_amplitudes)
      return
    end if

    this%outgoing_momenta(:, 1:2) = -momenta(:, 1:2)
    this%outgoing_momenta(:, 3:this%total_gluons) = &
         momenta(:, 3:this%total_gluons)

    ! Parke--Taylor gives every MHV or anti-MHV partial amplitude up to one
    ! event-wide phase, which cancels from the colour quadratic form.  When
    ! the complex vector is requested, evaluate only the canonical BG order
    ! and use it to restore the recursive external-wavefunction phase before
    ! filling the remaining DDM orders analytically.
    if (analytic_mhv_enabled .and. (outgoing_positive == 2 .or. &
         outgoing_positive == this%total_gluons-2)) then
      if (present(adjoint_amplitudes)) then
        do gluon = 1, this%total_gluons
          call external_massless_vector(this%outgoing_momenta(:, gluon), &
               helicities(gluon), this%external_wavefunctions(:, gluon))
        end do
      end if
      call evaluate_mhv_orders(this, helicities, coupling, outgoing_positive, &
           present(adjoint_amplitudes))
      call this%colour%contract_reusing_workspace( &
           this%amplitudes, matrix_element_squared)
      if (average_colours) matrix_element_squared = &
           matrix_element_squared/initial_colour_states
      if (present(adjoint_amplitudes)) &
           call copy_adjoint_amplitudes(this, adjoint_amplitudes)
      return
    end if

    do gluon = 1, this%total_gluons
      call external_massless_vector(this%outgoing_momenta(:, gluon), &
           helicities(gluon), this%external_wavefunctions(:, gluon))
    end do

    if (this%current_gluons <= maximum_shared_current_degree) then
      call ensure_current_dag(this)
      call this%current_dag%evaluate(this%outgoing_momenta(:, 1:this%current_gluons), &
           this%external_wavefunctions(:, 1:this%current_gluons), &
           this%external_wavefunctions(:, this%total_gluons), coupling, &
           this%amplitudes)
    else
      do order = 1, size(this%amplitudes)
        call evaluate_colour_order(this, order, coupling)
      end do
    end if
    call this%colour%contract_reusing_workspace( &
         this%amplitudes, matrix_element_squared)
    if (average_colours) &
         matrix_element_squared = matrix_element_squared/initial_colour_states

    if (present(adjoint_amplitudes)) then
      call copy_adjoint_amplitudes(this, adjoint_amplitudes)
    end if
  end subroutine evaluate_ampligluon_adjoint

  subroutine ensure_current_dag(this)
    class(ampligluon_adjoint_t), intent(inout) :: this

    if (allocated(this%current_dag)) return
    allocate(this%current_dag)
    call this%current_dag%initialize(this%current_gluons, &
         size(this%amplitudes))
  end subroutine ensure_current_dag

  subroutine copy_adjoint_amplitudes(this, adjoint_amplitudes)
    class(ampligluon_adjoint_t), intent(in) :: this
    complex(dp), allocatable, intent(inout) :: adjoint_amplitudes(:)

    if (allocated(adjoint_amplitudes)) then
      if (size(adjoint_amplitudes) /= size(this%amplitudes)) &
           deallocate(adjoint_amplitudes)
    end if
    if (.not. allocated(adjoint_amplitudes)) &
         allocate(adjoint_amplitudes(size(this%amplitudes)))
    adjoint_amplitudes = this%amplitudes
  end subroutine copy_adjoint_amplitudes

  subroutine evaluate_mhv_orders(this, helicities, coupling, outgoing_positive, &
       preserve_recursive_phase)
    class(ampligluon_adjoint_t), intent(inout) :: this
    integer, intent(in) :: helicities(:), outgoing_positive
    real(dp), intent(in) :: coupling
    logical, intent(in) :: preserve_recursive_phase

    complex(dp) :: angle_spinor(2, this%total_gluons)
    complex(dp) :: amplitude_normalization, bracket, canonical_amplitude
    complex(dp) :: inverse_denominator, numerator
    complex(dp) :: inverse_spinor_product(this%total_gluons, this%total_gluons)
    complex(dp) :: square_spinor(2, this%total_gluons)
    integer :: first, gluon, order, outgoing_helicity, pair(2), pair_count
    integer :: second
    logical :: use_square_brackets

    use_square_brackets = outgoing_positive == 2
    pair_count = 0
    do gluon = 1, this%total_gluons
      call massless_spinors(this%outgoing_momenta(:, gluon), &
           angle_spinor(:, gluon), square_spinor(:, gluon))
      outgoing_helicity = helicities(gluon)
      if (gluon <= 2) outgoing_helicity = -outgoing_helicity
      if ((use_square_brackets .and. outgoing_helicity > 0) .or. &
           (.not. use_square_brackets .and. outgoing_helicity < 0)) then
        pair_count = pair_count+1
        if (pair_count <= size(pair)) pair(pair_count) = gluon
      end if
    end do
    if (pair_count /= 2) call fail('wrong MHV distinguished-helicity count')

    inverse_spinor_product = cmplx(0.0_dp, 0.0_dp, dp)
    do first = 1, this%total_gluons-1
      do second = first+1, this%total_gluons
        if (use_square_brackets) then
          bracket = spinor_bracket( &
               square_spinor(:, first), square_spinor(:, second))
        else
          bracket = spinor_bracket( &
               angle_spinor(:, first), angle_spinor(:, second))
        end if
        if (abs(bracket) <= tiny(1.0_dp)) &
             call fail('singular spinor bracket in MHV evaluation')
        inverse_spinor_product(first, second) = 1.0_dp/bracket
        inverse_spinor_product(second, first) = &
             -inverse_spinor_product(first, second)
      end do
    end do

    if (preserve_recursive_phase) then
      call evaluate_colour_order(this, 1, coupling)
      canonical_amplitude = this%amplitudes(1)
      inverse_denominator = inverse_cyclic_denominator( &
           inverse_spinor_product, this%permutations(:, 1), &
           this%total_gluons)
      amplitude_normalization = canonical_amplitude/inverse_denominator
      do order = 2, size(this%amplitudes)
        this%amplitudes(order) = amplitude_normalization* &
             inverse_cyclic_denominator(inverse_spinor_product, &
             this%permutations(:, order), this%total_gluons)
      end do
    else
      numerator = coupling**(this%total_gluons-2)/ &
           inverse_spinor_product(pair(1), pair(2))**4
      do order = 1, size(this%amplitudes)
        this%amplitudes(order) = numerator*inverse_cyclic_denominator( &
             inverse_spinor_product, this%permutations(:, order), &
             this%total_gluons)
      end do
    end if
  end subroutine evaluate_mhv_orders

  pure complex(dp) function inverse_cyclic_denominator(inverse_products, &
       permutation, terminal) result(value)
    complex(dp), intent(in) :: inverse_products(:, :)
    integer, intent(in) :: permutation(:), terminal

    integer :: position

    value = inverse_products(terminal, permutation(1))* &
         inverse_products(permutation(size(permutation)), terminal)
    do position = 1, size(permutation)-1
      value = value*inverse_products( &
           permutation(position), permutation(position+1))
    end do
  end function inverse_cyclic_denominator

  subroutine massless_spinors(momentum, angle_spinor, square_spinor)
    real(dp), intent(in) :: momentum(0:3)
    complex(dp), intent(out) :: angle_spinor(2), square_spinor(2)

    complex(dp) :: root
    real(dp) :: minus_component, plus_component

    plus_component = momentum(0)+momentum(3)
    minus_component = momentum(0)-momentum(3)
    if (abs(plus_component) >= abs(minus_component)) then
      root = sqrt(cmplx(plus_component, 0.0_dp, dp))
      if (abs(root) <= tiny(1.0_dp)) &
           call fail('cannot factorize a zero massless momentum')
      angle_spinor(1) = root
      angle_spinor(2) = cmplx(momentum(1), momentum(2), dp)/root
      square_spinor(1) = root
      square_spinor(2) = cmplx(momentum(1), -momentum(2), dp)/root
    else
      root = sqrt(cmplx(minus_component, 0.0_dp, dp))
      if (abs(root) <= tiny(1.0_dp)) &
           call fail('cannot factorize a zero massless momentum')
      angle_spinor(1) = cmplx(momentum(1), -momentum(2), dp)/root
      angle_spinor(2) = root
      square_spinor(1) = cmplx(momentum(1), momentum(2), dp)/root
      square_spinor(2) = root
    end if
  end subroutine massless_spinors

  pure complex(dp) function spinor_bracket(first, second) result(bracket)
    complex(dp), intent(in) :: first(2), second(2)

    bracket = first(1)*second(2)-first(2)*second(1)
  end function spinor_bracket

  subroutine evaluate_colour_order(this, order, coupling)
    class(ampligluon_adjoint_t), intent(inout) :: this
    integer, intent(in) :: order
    real(dp), intent(in) :: coupling

    complex(dp) :: tensor_vertex(6), vector_vertex(4)
    integer :: external_gluon, first, last, left_last, length, position
    logical :: build_tensor

    this%interval_momenta = 0.0_dp
    this%gluon_current = (0.0_dp, 0.0_dp)
    this%tensor_current = (0.0_dp, 0.0_dp)
    do position = 1, this%current_gluons
      external_gluon = this%permutations(position, order)
      this%interval_momenta(:, position, position) = &
           this%outgoing_momenta(:, external_gluon)
      this%gluon_current(:, position, position) = &
           this%external_wavefunctions(:, external_gluon)
    end do

    do length = 2, this%current_gluons
      build_tensor = length < this%current_gluons
      do first = 1, this%current_gluons-length+1
        last = first+length-1
        this%interval_momenta(:, first, last) = &
             this%interval_momenta(:, first, last-1)+ &
             this%interval_momenta(:, last, last)
        do left_last = first, last-1
          call three_gluon( &
               this%gluon_current(:, first, left_last), &
               this%interval_momenta(:, first, left_last), &
               this%gluon_current(:, left_last+1, last), &
               this%interval_momenta(:, left_last+1, last), vector_vertex)
          this%gluon_current(:, first, last) = &
               this%gluon_current(:, first, last)+coupling*vector_vertex

          if (left_last > first) then
            call aux_tensor_gluon_to_gluon( &
                 this%tensor_current(:, first, left_last), &
                 this%gluon_current(:, left_last+1, last), vector_vertex)
            this%gluon_current(:, first, last) = &
                 this%gluon_current(:, first, last)+coupling*vector_vertex
          end if
          if (left_last+1 < last) then
            call gluon_aux_tensor_to_gluon( &
                 this%gluon_current(:, first, left_last), &
                 this%tensor_current(:, left_last+1, last), vector_vertex)
            this%gluon_current(:, first, last) = &
                 this%gluon_current(:, first, last)+coupling*vector_vertex
          end if

          if (build_tensor) then
            call two_gluon_to_aux_tensor( &
                 this%gluon_current(:, first, left_last), &
                 this%gluon_current(:, left_last+1, last), tensor_vertex)
            this%tensor_current(:, first, last) = &
                 this%tensor_current(:, first, last)+coupling*tensor_vertex
          end if
        end do
        if (length < this%current_gluons) &
             call massless_vector_propagator( &
                  this%gluon_current(:, first, last), &
                  this%interval_momenta(:, first, last))
      end do
    end do

    this%amplitudes(order) = terminal_vector_contract( &
         this%gluon_current(:, 1, this%current_gluons), &
         this%external_wavefunctions(:, this%total_gluons))
  end subroutine evaluate_colour_order

  subroutine validate_event(this, momenta, helicities)
    class(ampligluon_adjoint_t), intent(in) :: this
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

  integer function number_of_basis_amplitudes(this) result(number)
    class(ampligluon_adjoint_t), intent(in) :: this

    if (.not. this%initialized) &
         call fail('amplitude object is not initialized')
    number = this%colour%number_of_basis_elements()
  end function number_of_basis_amplitudes

  integer function number_of_final_gluons(this) result(number)
    class(ampligluon_adjoint_t), intent(in) :: this

    if (.not. this%initialized) &
         call fail('amplitude object is not initialized')
    number = this%total_gluons-2
  end function number_of_final_gluons

end module ampligluon_adjoint
