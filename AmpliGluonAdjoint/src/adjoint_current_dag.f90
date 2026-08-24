module adjoint_current_dag
  use ampligluon_common, only: dp, fail
  use trace_current_dag, only: trace_current_dag_t
  implicit none
  private

  type, public :: adjoint_current_dag_t
    private
    logical :: initialized = .false.
    integer :: degree_value = 0
    integer :: number_of_orders_value = 0
    type(trace_current_dag_t) :: trace_dag
  contains
    procedure, public :: evaluate => evaluate_adjoint_current_dag
    procedure, public :: initialize => initialize_adjoint_current_dag
    procedure, public :: number_of_orders
  end type adjoint_current_dag_t

contains

  subroutine initialize_adjoint_current_dag(this, degree, number_of_orders)
    class(adjoint_current_dag_t), intent(inout) :: this
    integer, intent(in) :: degree, number_of_orders

    integer :: expected_number_of_orders

    if (degree < 3) &
         call fail('adjoint current DAG degree must be at least three')
    if (number_of_orders < 1) &
         call fail('adjoint current DAG needs positive order count')

    this%initialized = .false.
    this%degree_value = degree
    this%number_of_orders_value = number_of_orders

    call this%trace_dag%initialize(degree, fixed_first=.true.)
    expected_number_of_orders = factorial_32(degree-1)
    if (number_of_orders /= expected_number_of_orders) then
      call fail(&
           'adjoint current DAG order count does not match fixed-first factorial')
    end if
    if (this%trace_dag%number_of_orders() /= number_of_orders) &
         call fail('adjoint and fixed-first current order counts disagree')
    this%initialized = .true.
  end subroutine initialize_adjoint_current_dag

  subroutine evaluate_adjoint_current_dag(this, outgoing_momenta, &
       external_wavefunctions, terminal_wavefunction, coupling, amplitudes)
    class(adjoint_current_dag_t), intent(inout) :: this
    real(dp), intent(in) :: outgoing_momenta(0:, :)
    complex(dp), intent(in) :: external_wavefunctions(:, :)
    complex(dp), intent(in) :: terminal_wavefunction(:)
    real(dp), intent(in) :: coupling
    complex(dp), intent(out) :: amplitudes(:)

    if (.not. this%initialized) call fail('adjoint current DAG is not initialized')
    if (size(outgoing_momenta, 1) /= 4 .or. &
         size(outgoing_momenta, 2) /= this%degree_value) then
      call fail('wrong momentum array for adjoint current DAG')
    end if
    if (size(external_wavefunctions, 1) /= 4 .or. &
         size(external_wavefunctions, 2) /= this%degree_value) then
      call fail('wrong wavefunction array for adjoint current DAG')
    end if
    if (size(terminal_wavefunction) /= 4) then
      call fail('wrong terminal wavefunction for adjoint current DAG')
    end if
    if (size(amplitudes) /= this%number_of_orders_value) then
      call fail('wrong amplitude array for adjoint current DAG')
    end if

    call this%trace_dag%evaluate(outgoing_momenta, external_wavefunctions, &
         terminal_wavefunction, coupling, amplitudes)
  end subroutine evaluate_adjoint_current_dag

  integer function number_of_orders(this) result(number)
    class(adjoint_current_dag_t), intent(in) :: this

    if (.not. this%initialized) call fail('adjoint current DAG is not initialized')
    number = this%number_of_orders_value
  end function number_of_orders

  integer function factorial_32(number) result(value)
    integer, intent(in) :: number

    integer :: factor

    value = 1
    if (number < 0) call fail('adjoint factorial requires non-negative input')
    do factor = 2, number
      if (value > huge(value)/factor) call fail('adjoint factorial overflows integer')
      value = value*factor
    end do
  end function factorial_32

end module adjoint_current_dag
