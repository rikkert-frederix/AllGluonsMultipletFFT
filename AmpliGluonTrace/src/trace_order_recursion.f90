module trace_order_recursion
  use ampligluon_common, only: dp, fail, i64
  use gluon_kinematics, only: aux_tensor_gluon_to_gluon, &
       gluon_aux_tensor_to_gluon, massless_vector_propagator, &
       terminal_vector_contract, three_gluon, two_gluon_to_aux_tensor
  implicit none
  private

  public :: evaluate_canonical_trace_order

  type, public :: trace_order_recursion_t
    private
    logical :: initialized = .false.
    integer :: degree_value = 0
    integer :: number_of_orders_value = 0
    integer, allocatable :: permutations(:, :)
    integer, allocatable :: reflection_partner(:)
    real(dp), allocatable :: interval_momenta(:, :, :)
    complex(dp), allocatable :: gluon_current(:, :, :)
    complex(dp), allocatable :: tensor_current(:, :, :)
  contains
    procedure, public :: evaluate => evaluate_trace_order_recursion
    procedure, public :: initialize => initialize_trace_order_recursion
    procedure, public :: number_of_orders
  end type trace_order_recursion_t

contains

  subroutine initialize_trace_order_recursion(this, degree, number_of_orders)
    class(trace_order_recursion_t), intent(inout) :: this
    integer, intent(in) :: degree, number_of_orders

    integer :: order, partner, position

    if (degree < 3) call fail('trace order degree must be at least three')
    if (factorial_i64(degree) /= int(number_of_orders, i64)) &
         call fail('trace order recursion has the wrong number of orders')
    this%initialized = .false.
    this%degree_value = degree
    this%number_of_orders_value = number_of_orders
    if (allocated(this%permutations)) deallocate(this%permutations)
    if (allocated(this%reflection_partner)) &
         deallocate(this%reflection_partner)
    if (allocated(this%interval_momenta)) deallocate(this%interval_momenta)
    if (allocated(this%gluon_current)) deallocate(this%gluon_current)
    if (allocated(this%tensor_current)) deallocate(this%tensor_current)

    allocate(this%permutations(degree, number_of_orders))
    this%permutations(:, 1) = [(position, position=1, degree)]
    do order = 2, number_of_orders
      this%permutations(:, order) = this%permutations(:, order-1)
      call next_permutation(this%permutations(:, order))
    end do
    allocate(this%reflection_partner(number_of_orders))
    do order = 1, number_of_orders
      partner = lexicographic_rank( &
           this%permutations(degree:1:-1, order))
      if (partner < 1 .or. partner > number_of_orders) &
           call fail('reflected trace order has an invalid rank')
      this%reflection_partner(order) = partner
    end do
    do order = 1, number_of_orders
      partner = this%reflection_partner(order)
      if (partner == order .or. &
           this%reflection_partner(partner) /= order) &
           call fail('trace-order reflection is not a free involution')
    end do

    allocate(this%interval_momenta(0:3, degree, degree))
    allocate(this%gluon_current(4, degree, degree))
    allocate(this%tensor_current(6, degree, degree))
    this%initialized = .true.
  end subroutine initialize_trace_order_recursion

  subroutine evaluate_trace_order_recursion(this, outgoing_momenta, &
       external_wavefunctions, terminal_wavefunction, coupling, amplitudes)
    class(trace_order_recursion_t), intent(inout) :: this
    real(dp), intent(in) :: outgoing_momenta(0:, :)
    complex(dp), intent(in) :: external_wavefunctions(:, :)
    complex(dp), intent(in) :: terminal_wavefunction(:)
    real(dp), intent(in) :: coupling
    complex(dp), intent(out) :: amplitudes(:)

    real(dp) :: reflection_sign
    integer :: order, partner

    if (.not. this%initialized) &
         call fail('trace order recursion is not initialized')
    if (size(outgoing_momenta, 1) /= 4 .or. &
         size(outgoing_momenta, 2) /= this%degree_value) &
         call fail('wrong momentum array for trace order recursion')
    if (size(external_wavefunctions, 1) /= 4 .or. &
         size(external_wavefunctions, 2) /= this%degree_value) &
         call fail('wrong wavefunction array for trace order recursion')
    if (size(terminal_wavefunction) /= 4) &
         call fail('wrong terminal wavefunction for trace order recursion')
    if (size(amplitudes) /= this%number_of_orders_value) &
         call fail('wrong amplitude array for trace order recursion')

    reflection_sign = merge(1.0_dp, -1.0_dp, &
         modulo(this%degree_value+1, 2) == 0)
    do order = 1, this%number_of_orders_value
      partner = this%reflection_partner(order)
      if (order > partner) cycle
      call evaluate_one_order(this, order, outgoing_momenta, &
           external_wavefunctions, terminal_wavefunction, coupling, &
           amplitudes(order))
      amplitudes(partner) = reflection_sign*amplitudes(order)
    end do
  end subroutine evaluate_trace_order_recursion

  subroutine evaluate_one_order(this, order, outgoing_momenta, &
       external_wavefunctions, terminal_wavefunction, coupling, amplitude)
    class(trace_order_recursion_t), intent(inout) :: this
    integer, intent(in) :: order
    real(dp), intent(in) :: outgoing_momenta(0:, :)
    complex(dp), intent(in) :: external_wavefunctions(:, :)
    complex(dp), intent(in) :: terminal_wavefunction(:)
    real(dp), intent(in) :: coupling
    complex(dp), intent(out) :: amplitude

    call evaluate_trace_order_with_workspace(this%permutations(:, order), &
         outgoing_momenta, external_wavefunctions, terminal_wavefunction, &
         coupling, this%interval_momenta, this%gluon_current, &
         this%tensor_current, amplitude)
  end subroutine evaluate_one_order

  subroutine evaluate_canonical_trace_order(outgoing_momenta, &
       external_wavefunctions, terminal_wavefunction, coupling, amplitude)
    real(dp), intent(in) :: outgoing_momenta(0:, :)
    complex(dp), intent(in) :: external_wavefunctions(:, :)
    complex(dp), intent(in) :: terminal_wavefunction(:)
    real(dp), intent(in) :: coupling
    complex(dp), intent(out) :: amplitude

    integer :: degree, position
    integer :: permutation(size(outgoing_momenta, 2))
    real(dp) :: interval_momenta(0:3, size(outgoing_momenta, 2), &
         size(outgoing_momenta, 2))
    complex(dp) :: gluon_current(4, size(outgoing_momenta, 2), &
         size(outgoing_momenta, 2))
    complex(dp) :: tensor_current(6, size(outgoing_momenta, 2), &
         size(outgoing_momenta, 2))

    degree = size(outgoing_momenta, 2)
    if (degree < 3 .or. size(outgoing_momenta, 1) /= 4) &
         call fail('wrong momentum array for canonical trace order')
    if (size(external_wavefunctions, 1) /= 4 .or. &
         size(external_wavefunctions, 2) /= degree) &
         call fail('wrong wavefunction array for canonical trace order')
    if (size(terminal_wavefunction) /= 4) &
         call fail('wrong terminal wavefunction for canonical trace order')
    permutation = [(position, position=1, degree)]
    call evaluate_trace_order_with_workspace(permutation, outgoing_momenta, &
         external_wavefunctions, terminal_wavefunction, coupling, &
         interval_momenta, gluon_current, tensor_current, amplitude)
  end subroutine evaluate_canonical_trace_order

  subroutine evaluate_trace_order_with_workspace(permutation, &
       outgoing_momenta, external_wavefunctions, terminal_wavefunction, &
       coupling, interval_momenta, gluon_current, tensor_current, amplitude)
    integer, intent(in) :: permutation(:)
    real(dp), intent(in) :: outgoing_momenta(0:, :)
    complex(dp), intent(in) :: external_wavefunctions(:, :)
    complex(dp), intent(in) :: terminal_wavefunction(:)
    real(dp), intent(in) :: coupling
    real(dp), intent(inout) :: interval_momenta(0:, :, :)
    complex(dp), intent(inout) :: gluon_current(:, :, :)
    complex(dp), intent(inout) :: tensor_current(:, :, :)
    complex(dp), intent(out) :: amplitude

    complex(dp) :: tensor_vertex(6), vector_vertex(4)
    integer :: degree, external_gluon, first, last, left_last, length, position

    degree = size(permutation)
    interval_momenta = 0.0_dp
    gluon_current = cmplx(0.0_dp, 0.0_dp, dp)
    tensor_current = cmplx(0.0_dp, 0.0_dp, dp)
    do position = 1, degree
      external_gluon = permutation(position)
      interval_momenta(:, position, position) = &
           outgoing_momenta(:, external_gluon)
      gluon_current(:, position, position) = &
           external_wavefunctions(:, external_gluon)
    end do

    do length = 2, degree
      do first = 1, degree-length+1
        last = first+length-1
        interval_momenta(:, first, last) = &
             interval_momenta(:, first, last-1)+ &
             interval_momenta(:, last, last)
        do left_last = first, last-1
          call three_gluon(gluon_current(:, first, left_last), &
               interval_momenta(:, first, left_last), &
               gluon_current(:, left_last+1, last), &
               interval_momenta(:, left_last+1, last), vector_vertex)
          gluon_current(:, first, last) = &
               gluon_current(:, first, last)+coupling*vector_vertex
          if (left_last > first) then
            call aux_tensor_gluon_to_gluon( &
                 tensor_current(:, first, left_last), &
                 gluon_current(:, left_last+1, last), vector_vertex)
            gluon_current(:, first, last) = &
                 gluon_current(:, first, last)+coupling*vector_vertex
          end if
          if (left_last+1 < last) then
            call gluon_aux_tensor_to_gluon( &
                 gluon_current(:, first, left_last), &
                 tensor_current(:, left_last+1, last), vector_vertex)
            gluon_current(:, first, last) = &
                 gluon_current(:, first, last)+coupling*vector_vertex
          end if
          if (length < degree) then
            call two_gluon_to_aux_tensor( &
                 gluon_current(:, first, left_last), &
                 gluon_current(:, left_last+1, last), tensor_vertex)
            tensor_current(:, first, last) = &
                 tensor_current(:, first, last)+coupling*tensor_vertex
          end if
        end do
        if (length < degree) call massless_vector_propagator( &
             gluon_current(:, first, last), &
             interval_momenta(:, first, last))
      end do
    end do
    amplitude = terminal_vector_contract( &
         gluon_current(:, 1, degree), terminal_wavefunction)
  end subroutine evaluate_trace_order_with_workspace

  integer function lexicographic_rank(permutation) result(rank)
    integer, intent(in) :: permutation(:)

    integer :: factor, left, right, smaller

    factor = 1
    do left = 2, size(permutation)-1
      if (factor > huge(factor)/left) &
           call fail('trace-order rank exceeds integer capacity')
      factor = factor*left
    end do
    rank = 1
    do left = 1, size(permutation)-1
      smaller = 0
      do right = left+1, size(permutation)
        if (permutation(right) < permutation(left)) smaller = smaller+1
      end do
      rank = rank+smaller*factor
      factor = factor/(size(permutation)-left)
    end do
  end function lexicographic_rank

  subroutine next_permutation(permutation)
    integer, intent(inout) :: permutation(:)

    integer :: left, right, temporary

    left = size(permutation)-1
    do while (left >= 1)
      if (permutation(left) < permutation(left+1)) exit
      left = left-1
    end do
    if (left < 1) call fail('trace-order permutation wrapped unexpectedly')
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
           call fail('trace-order factorial overflows 64 bits')
      value = value*int(factor, i64)
    end do
  end function factorial_i64

  integer function number_of_orders(this) result(number)
    class(trace_order_recursion_t), intent(in) :: this

    if (.not. this%initialized) &
         call fail('trace order recursion is not initialized')
    number = this%number_of_orders_value
  end function number_of_orders

end module trace_order_recursion
