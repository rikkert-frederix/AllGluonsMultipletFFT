module trace_colour_kernel
  use ampligluon_common, only: dp, fail, i64
  implicit none
  private

  integer(i64), parameter :: fundamental_colours = 3_i64

  public :: build_u3_trace_kernel

contains

  subroutine build_u3_trace_kernel(total_gluons, kernel)
    integer, intent(in) :: total_gluons
    real(dp), allocatable, intent(out) :: kernel(:)

    integer, allocatable :: commutator(:), inverse(:), permutation(:)
    integer, allocatable :: visit_mark(:)
    integer(i64) :: kernel_value, number_of_orders_64
    integer :: cycle_count, label, order, position

    if (total_gluons < 4) &
         call fail('trace colour kernel needs at least four total gluons')
    number_of_orders_64 = factorial_i64(total_gluons-1)
    if (number_of_orders_64 > int(huge(0), i64)) &
         call fail('trace colour kernel exceeds array-index capacity')

    if (allocated(kernel)) deallocate(kernel)
    allocate(kernel(int(number_of_orders_64)))
    allocate(permutation(total_gluons), inverse(total_gluons))
    allocate(commutator(total_gluons), visit_mark(total_gluons))
    visit_mark = 0
    permutation = [(position, position=1, total_gluons)]

    do order = 1, size(kernel)
      do position = 1, total_gluons
        inverse(permutation(position)) = position
      end do

      ! Embed h in S_N by fixing N.  For c=(1 2 ... N), the U(3) trace
      ! kernel is 3 to the number of cycles of c^{-1} h c h^{-1}.
      do label = 1, total_gluons
        position = inverse(label)
        position = modulo(position, total_gluons)+1
        position = permutation(position)
        commutator(label) = modulo(position-2, total_gluons)+1
      end do
      cycle_count = number_of_cycles(commutator, visit_mark, order)
      kernel_value = integer_power(fundamental_colours, cycle_count)
      kernel(order) = real(kernel_value, dp)

      if (order < size(kernel)) &
           call next_permutation(permutation(1:total_gluons-1))
    end do
  end subroutine build_u3_trace_kernel

  integer function number_of_cycles(permutation, visit_mark, mark) &
       result(number)
    integer, intent(in) :: permutation(:)
    integer, intent(inout) :: visit_mark(:)
    integer, intent(in) :: mark

    integer :: item, position

    if (size(visit_mark) /= size(permutation) .or. mark < 1) &
         call fail('invalid trace-kernel cycle-count workspace')
    number = 0
    do item = 1, size(permutation)
      if (visit_mark(item) == mark) cycle
      number = number+1
      position = item
      do while (visit_mark(position) /= mark)
        visit_mark(position) = mark
        position = permutation(position)
        if (position < 1 .or. position > size(permutation)) &
             call fail('invalid permutation in trace colour kernel')
      end do
    end do
  end function number_of_cycles

  integer(i64) function factorial_i64(number) result(value)
    integer, intent(in) :: number

    integer :: factor

    if (number < 0) call fail('negative trace-kernel factorial')
    value = 1_i64
    do factor = 2, number
      if (value > huge(value)/int(factor, i64)) &
           call fail('trace-kernel factorial overflows a 64-bit integer')
      value = value*int(factor, i64)
    end do
  end function factorial_i64

  integer(i64) function integer_power(base, exponent) result(value)
    integer(i64), intent(in) :: base
    integer, intent(in) :: exponent

    integer :: power

    if (base < 1_i64 .or. exponent < 0) &
         call fail('invalid integer power in trace colour kernel')
    value = 1_i64
    do power = 1, exponent
      if (value > huge(value)/base) &
           call fail('trace colour kernel overflows a 64-bit integer')
      value = value*base
    end do
  end function integer_power

  subroutine next_permutation(permutation)
    integer, intent(inout) :: permutation(:)

    integer :: left, right, temporary

    left = size(permutation)-1
    do while (left >= 1)
      if (permutation(left) < permutation(left+1)) exit
      left = left-1
    end do
    if (left < 1) &
         call fail('trace-kernel permutation sequence wrapped unexpectedly')
    right = size(permutation)
    do while (permutation(right) < permutation(left))
      right = right-1
    end do
    temporary = permutation(left)
    permutation(left) = permutation(right)
    permutation(right) = temporary
    permutation(left+1:) = permutation(size(permutation):left+1:-1)
  end subroutine next_permutation

end module trace_colour_kernel
