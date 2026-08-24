module trace_mhv
  use ampligluon_common, only: dp, fail, i64
  use trace_order_recursion, only: evaluate_canonical_trace_order
  implicit none
  private

  type, public :: trace_mhv_t
    private
    logical :: initialized = .false.
    integer :: total_gluons = 0
    integer :: number_of_orders_value = 0
    complex(dp), allocatable :: inverse_bracket(:, :)
    integer, allocatable :: permutation(:)
    integer, allocatable :: reflection_partner(:)
  contains
    procedure, public :: evaluate => evaluate_trace_mhv
    procedure, public :: initialize => initialize_trace_mhv
  end type trace_mhv_t

contains

  subroutine initialize_trace_mhv(this, total_gluons, number_of_orders)
    class(trace_mhv_t), intent(inout) :: this
    integer, intent(in) :: total_gluons, number_of_orders

    integer(i64) :: expected_orders
    integer(i64) :: factorials(0:total_gluons-1)
    integer :: order, position

    if (total_gluons < 4) &
         call fail('trace MHV evaluation needs at least four gluons')
    expected_orders = factorial_i64(total_gluons-1)
    if (expected_orders /= int(number_of_orders, i64)) &
         call fail('trace MHV evaluator has the wrong number of orders')
    this%initialized = .false.
    this%total_gluons = total_gluons
    this%number_of_orders_value = number_of_orders
    if (allocated(this%inverse_bracket)) deallocate(this%inverse_bracket)
    if (allocated(this%permutation)) deallocate(this%permutation)
    if (allocated(this%reflection_partner)) &
         deallocate(this%reflection_partner)
    allocate(this%inverse_bracket(total_gluons, total_gluons))
    allocate(this%permutation(total_gluons-1))
    allocate(this%reflection_partner(number_of_orders))
    factorials(0) = 1_i64
    do position = 1, total_gluons-1
      factorials(position) = int(position, i64)*factorials(position-1)
    end do
    this%permutation = [(position, position=1, total_gluons-1)]
    do order = 1, number_of_orders
      this%reflection_partner(order) = &
           reversed_permutation_rank(this%permutation, factorials)
      if (order < number_of_orders) call next_permutation(this%permutation)
    end do
    this%initialized = .true.
  end subroutine initialize_trace_mhv

  subroutine evaluate_trace_mhv(this, outgoing_momenta, &
       external_wavefunctions, physical_helicities, coupling, amplitudes)
    class(trace_mhv_t), intent(inout) :: this
    real(dp), intent(in) :: outgoing_momenta(0:, :)
    complex(dp), intent(in) :: external_wavefunctions(:, :)
    integer, intent(in) :: physical_helicities(:)
    real(dp), intent(in) :: coupling
    complex(dp), intent(out) :: amplitudes(:)

    complex(dp) :: amplitude_normalization, canonical_amplitude
    integer :: order, outgoing_positive, position
    real(dp) :: reflection_sign
    logical :: use_angle_brackets

    if (.not. this%initialized) &
         call fail('trace MHV evaluator is not initialized')
    if (size(outgoing_momenta, 1) /= 4 .or. &
         size(outgoing_momenta, 2) /= this%total_gluons) &
         call fail('wrong momentum array for trace MHV evaluation')
    if (size(external_wavefunctions, 1) /= 4 .or. &
         size(external_wavefunctions, 2) /= this%total_gluons) &
         call fail('wrong wavefunction array for trace MHV evaluation')
    if (size(physical_helicities) /= this%total_gluons) &
         call fail('wrong helicity array for trace MHV evaluation')
    if (size(amplitudes) /= this%number_of_orders_value) &
         call fail('wrong amplitude array for trace MHV evaluation')

    outgoing_positive = count(physical_helicities(1:2) < 0)+ &
         count(physical_helicities(3:) > 0)
    if (outgoing_positive /= 2 .and. &
         outgoing_positive /= this%total_gluons-2) &
      call fail('trace MHV evaluator received a non-MHV helicity sector')
    ! At four points both descriptions apply; square brackets are sufficient.
    use_angle_brackets = outgoing_positive == this%total_gluons-2 .and. &
         outgoing_positive /= 2

    call build_inverse_spinor_brackets(outgoing_momenta, use_angle_brackets, &
         this%inverse_bracket)
    this%permutation = [(position, position=1, this%total_gluons-1)]
    call evaluate_canonical_trace_order( &
         outgoing_momenta(:, :this%total_gluons-1), &
         external_wavefunctions(:, :this%total_gluons-1), &
         external_wavefunctions(:, this%total_gluons), coupling, &
         canonical_amplitude)

    amplitude_normalization = canonical_amplitude/ &
         inverse_cyclic_denominator(this%inverse_bracket, this%permutation, &
         this%total_gluons)
    amplitudes(1) = canonical_amplitude
    reflection_sign = 1.0_dp
    if (modulo(this%total_gluons, 2) /= 0) reflection_sign = -1.0_dp
    do order = 2, this%number_of_orders_value
      call next_permutation(this%permutation)
      if (order > this%reflection_partner(order)) then
        amplitudes(order) = reflection_sign* &
             amplitudes(this%reflection_partner(order))
      else
        amplitudes(order) = amplitude_normalization* &
             inverse_cyclic_denominator(this%inverse_bracket, &
             this%permutation, this%total_gluons)
      end if
    end do
  end subroutine evaluate_trace_mhv

  integer function reversed_permutation_rank(permutation, factorials) &
       result(rank)
    integer, intent(in) :: permutation(:)
    integer(i64), intent(in) :: factorials(0:)

    integer(i64) :: lower_values, rank_64, used_values
    integer :: position, smaller, value

    if (size(factorials) <= size(permutation)) &
         call fail('trace MHV rank table is too short')
    rank_64 = 1_i64
    used_values = 0_i64
    do position = 1, size(permutation)
      value = permutation(size(permutation)-position+1)
      if (value == 1) then
        smaller = 0
      else
        lower_values = shiftl(1_i64, value-1)-1_i64
        smaller = popcnt(iand(lower_values, not(used_values)))
      end if
      rank_64 = rank_64+int(smaller, i64)* &
           factorials(size(permutation)-position)
      used_values = ibset(used_values, value-1)
    end do
    if (rank_64 > int(huge(rank), i64)) &
         call fail('trace MHV reflected rank exceeds array-index capacity')
    rank = int(rank_64)
  end function reversed_permutation_rank

  subroutine build_inverse_spinor_brackets(momenta, use_angle_brackets, &
       inverse_brackets)
    real(dp), intent(in) :: momenta(0:, :)
    logical, intent(in) :: use_angle_brackets
    complex(dp), intent(out) :: inverse_brackets(:, :)

    complex(dp) :: bracket, root
    complex(dp) :: spinor(2, size(momenta, 2))
    integer :: first, leg, second

    if (size(inverse_brackets, 1) /= size(momenta, 2) .or. &
         size(inverse_brackets, 2) /= size(momenta, 2)) &
         call fail('wrong spinor-bracket workspace shape')
    do leg = 1, size(momenta, 2)
      if (abs(momenta(0, leg)+momenta(3, leg)) > &
           abs(momenta(0, leg)-momenta(3, leg))) then
        root = sqrt(cmplx(momenta(0, leg)+momenta(3, leg), 0.0_dp, dp))
        spinor(1, leg) = root
        if (use_angle_brackets) then
          spinor(2, leg) = cmplx(momenta(1, leg), &
               momenta(2, leg), dp)/root
        else
          spinor(2, leg) = cmplx(momenta(1, leg), &
               -momenta(2, leg), dp)/root
        end if
      else
        root = sqrt(cmplx(momenta(0, leg)-momenta(3, leg), 0.0_dp, dp))
        spinor(2, leg) = root
        if (use_angle_brackets) then
          spinor(1, leg) = cmplx(momenta(1, leg), &
               -momenta(2, leg), dp)/root
        else
          spinor(1, leg) = cmplx(momenta(1, leg), &
               momenta(2, leg), dp)/root
        end if
      end if
    end do

    inverse_brackets = cmplx(0.0_dp, 0.0_dp, dp)
    do first = 1, size(momenta, 2)-1
      do second = first+1, size(momenta, 2)
        bracket = spinor(1, first)*spinor(2, second)- &
             spinor(2, first)*spinor(1, second)
        if (abs(bracket) <= tiny(1.0_dp)) &
             call fail('singular spinor bracket in trace MHV evaluation')
        inverse_brackets(first, second) = 1.0_dp/bracket
        inverse_brackets(second, first) = -inverse_brackets(first, second)
      end do
    end do
  end subroutine build_inverse_spinor_brackets

  complex(dp) function inverse_cyclic_denominator(inverse_brackets, &
       permutation, terminal) result(inverse_denominator)
    complex(dp), intent(in) :: inverse_brackets(:, :)
    integer, intent(in) :: permutation(:), terminal

    integer :: position

    inverse_denominator = inverse_brackets(terminal, permutation(1))* &
         inverse_brackets(permutation(size(permutation)), terminal)
    do position = 1, size(permutation)-1
      inverse_denominator = inverse_denominator* &
           inverse_brackets(permutation(position), permutation(position+1))
    end do
  end function inverse_cyclic_denominator

  subroutine next_permutation(permutation)
    integer, intent(inout) :: permutation(:)

    integer :: left, right, temporary

    left = size(permutation)-1
    do while (left >= 1)
      if (permutation(left) < permutation(left+1)) exit
      left = left-1
    end do
    if (left < 1) &
         call fail('trace MHV permutation sequence wrapped unexpectedly')
    right = size(permutation)
    do while (permutation(right) < permutation(left))
      right = right-1
    end do
    temporary = permutation(left)
    permutation(left) = permutation(right)
    permutation(right) = temporary
    permutation(left+1:) = permutation(size(permutation):left+1:-1)
  end subroutine next_permutation

  integer(i64) function factorial_i64(number) result(value)
    integer, intent(in) :: number

    integer :: factor

    value = 1_i64
    do factor = 2, number
      if (value > huge(value)/int(factor, i64)) &
           call fail('trace MHV factorial overflows 64 bits')
      value = value*int(factor, i64)
    end do
  end function factorial_i64

end module trace_mhv
