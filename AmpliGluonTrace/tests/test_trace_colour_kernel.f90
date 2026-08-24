program test_trace_colour_kernel
  use ampligluon_common, only: dp, fail, i64
  use trace_colour_kernel, only: build_u3_trace_kernel
  implicit none

  call check_kernel(4, [81_i64, 9_i64], [1, 5])
  call check_kernel(5, [243_i64, 27_i64, 3_i64], [1, 15, 8])
  call check_kernel(6, [729_i64, 81_i64, 9_i64], [1, 35, 84])
  call check_kernel(7, [2187_i64, 243_i64, 27_i64, 3_i64], &
                    [1, 70, 469, 180])
  call check_kernel(8, [6561_i64, 729_i64, 81_i64, 9_i64], &
                    [1, 126, 1869, 3044])
  write(*, '(a)') 'trace colour kernel regression: PASS'

contains

  subroutine check_kernel(total_gluons, values, multiplicities)
    integer, intent(in) :: total_gluons
    integer(i64), intent(in) :: values(:)
    integer, intent(in) :: multiplicities(:)

    real(dp), allocatable :: kernel(:)
    integer, allocatable :: inverse(:), permutation(:)
    integer :: factor, order, position

    call build_u3_trace_kernel(total_gluons, kernel)
    call require(size(values) == size(multiplicities), &
                 'kernel histogram has inconsistent dimensions')
    call require(size(kernel) == factorial(total_gluons-1), &
                 'kernel has the wrong group order')
    do factor = 1, size(values)
      call require(count(nint(kernel, i64) == values(factor)) == &
                   multiplicities(factor), &
                   'kernel value multiplicity changed')
    end do
    call require(sum(multiplicities) == size(kernel), &
                 'kernel histogram does not cover the group')

    allocate(permutation(total_gluons-1), inverse(total_gluons-1))
    permutation = [(position, position=1, total_gluons-1)]
    do order = 1, size(kernel)
      do position = 1, size(permutation)
        inverse(permutation(position)) = position
      end do
      call require(nint(kernel(order), i64) == &
                   nint(kernel(permutation_rank(inverse)), i64), &
                   'trace kernel is not invariant under inversion')
      if (order < size(kernel)) call next_permutation(permutation)
    end do
  end subroutine check_kernel

  integer function factorial(number) result(value)
    integer, intent(in) :: number
    integer :: factor

    value = 1
    do factor = 2, number
      if (value > huge(value)/factor) call fail('test factorial overflow')
      value = value*factor
    end do
  end function factorial

  integer function permutation_rank(permutation) result(rank)
    integer, intent(in) :: permutation(:)

    integer :: factor, left, right, smaller

    rank = 1
    factor = factorial(size(permutation)-1)
    do left = 1, size(permutation)-1
      smaller = 0
      do right = left+1, size(permutation)
        if (permutation(right) < permutation(left)) smaller = smaller+1
      end do
      rank = rank+smaller*factor
      factor = factor/(size(permutation)-left)
    end do
  end function permutation_rank

  subroutine next_permutation(permutation)
    integer, intent(inout) :: permutation(:)

    integer :: left, right, temporary

    left = size(permutation)-1
    do while (left >= 1)
      if (permutation(left) < permutation(left+1)) exit
      left = left-1
    end do
    if (left < 1) call fail('test permutation wrapped unexpectedly')
    right = size(permutation)
    do while (permutation(right) < permutation(left))
      right = right-1
    end do
    temporary = permutation(left)
    permutation(left) = permutation(right)
    permutation(right) = temporary
    permutation(left+1:) = permutation(size(permutation):left+1:-1)
  end subroutine next_permutation

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) call fail(message)
  end subroutine require

end program test_trace_colour_kernel
