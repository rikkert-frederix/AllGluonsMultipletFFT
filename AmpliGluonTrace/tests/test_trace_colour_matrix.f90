program test_trace_colour_matrix
  use ampligluon_common, only: dp, fail
  use trace_colour_kernel, only: build_u3_trace_kernel
  use trace_colour_matrix, only: trace_colour_matrix_t
  implicit none

  integer :: total_gluons

  do total_gluons = 4, 7
    call check_arbitrary_vector(total_gluons)
  end do
  call check_mode_reinitialization()
  write(*, '(a)') 'trace colour contraction regression: PASS'

contains

  subroutine check_arbitrary_vector(total_gluons)
    integer, intent(in) :: total_gluons

    type(trace_colour_matrix_t) :: direct_colour, fourier_colour
    complex(dp), allocatable :: amplitudes(:)
    integer, allocatable :: factorials(:), permutations(:, :)
    real(dp), allocatable :: kernel(:)
    complex(dp) :: direct
    real(dp) :: direct_result, fast, scale
    integer :: order, position, row

    call fourier_colour%initialize(total_gluons, permutations)
    call direct_colour%initialize(total_gluons, use_colour_fft=.false.)
    order = fourier_colour%number_of_orders()
    call require(order == factorial(total_gluons-1), &
         'colour contraction has the wrong group order')
    call require(direct_colour%number_of_orders() == order, &
         'direct colour contraction has the wrong group order')
    call require(fourier_colour%number_of_stored_kernel_values() == order, &
         'Fourier colour kernel storage is not linear in the group order')
    call require(direct_colour%number_of_stored_kernel_values() == order, &
         'direct colour kernel storage is not linear in the group order')
    allocate(amplitudes(order))
    do row = 1, order
      amplitudes(row) = cmplx( &
           sin(0.173_dp*real(row, dp))+0.003_dp*real(row, dp), &
           cos(0.119_dp*real(row, dp))-0.002_dp*real(row, dp), dp)
    end do
    call fourier_colour%contract(amplitudes, fast)
    call direct_colour%contract(amplitudes, direct_result)

    call build_u3_trace_kernel(total_gluons, kernel)
    allocate(factorials(0:total_gluons-1))
    factorials(0) = 1
    do position = 1, total_gluons-1
      factorials(position) = position*factorials(position-1)
    end do
    direct = direct_colour_sum(amplitudes, kernel, permutations, factorials)

    scale = max(1.0_dp, abs(real(direct, dp)))
    if (abs(aimag(direct)) > 2.0e-12_dp*scale .or. &
         abs(fast-real(direct, dp)) > 2.0e-10_dp*scale) then
      write(*, '(a,i0,3(1x,es24.16))') &
           'Fourier/direct mismatch for total gluons ', total_gluons, &
           fast, real(direct, dp), aimag(direct)
      call fail('trace colour Fourier contraction is incorrect')
    end if
    if (abs(direct_result-real(direct, dp)) > 2.0e-12_dp*scale) then
      write(*, '(a,i0,2(1x,es24.16))') &
           'normal/direct mismatch for total gluons ', total_gluons, &
           direct_result, real(direct, dp)
      call fail('normal trace colour contraction is incorrect')
    end if

    ! Reuse the same Trace object with an unrelated vector.  This exercises
    ! the persistent Fourier workspace on a second call.
    do row = 1, order
      amplitudes(row) = cmplx( &
           cos(0.151_dp*real(row, dp))-0.004_dp*real(row, dp), &
           sin(0.107_dp*real(row, dp))+0.006_dp*real(row, dp), dp)
    end do
    call fourier_colour%contract(amplitudes, fast)
    call direct_colour%contract(amplitudes, direct_result)
    direct = direct_colour_sum(amplitudes, kernel, permutations, factorials)
    scale = max(1.0_dp, abs(real(direct, dp)))
    if (abs(aimag(direct)) > 2.0e-12_dp*scale .or. &
         abs(fast-real(direct, dp)) > 2.0e-10_dp*scale) &
         call fail('reused trace colour workspace retained state')
    if (abs(direct_result-real(direct, dp)) > 2.0e-12_dp*scale) &
         call fail('reused normal trace colour contraction retained state')
  end subroutine check_arbitrary_vector

  subroutine check_mode_reinitialization()
    type(trace_colour_matrix_t) :: colour
    complex(dp) :: amplitudes(24)
    real(dp) :: direct_result, fourier_result, scale
    integer :: order

    do order = 1, size(amplitudes)
      amplitudes(order) = cmplx(sin(0.271_dp*real(order, dp)), &
           cos(0.193_dp*real(order, dp)), dp)
    end do

    call colour%initialize(5)
    call colour%contract(amplitudes, fourier_result)
    call colour%initialize(5, use_colour_fft=.false.)
    call colour%contract(amplitudes, direct_result)
    scale = max(1.0_dp, abs(fourier_result), abs(direct_result))
    call require(abs(fourier_result-direct_result) <= 2.0e-10_dp*scale, &
         'FFT-to-direct colour reinitialization changed the result')

    call colour%initialize(5)
    call colour%contract(amplitudes, fourier_result)
    call require(abs(fourier_result-direct_result) <= 2.0e-10_dp*scale, &
         'direct-to-FFT colour reinitialization changed the result')
  end subroutine check_mode_reinitialization

  function direct_colour_sum(amplitudes, kernel, permutations, factorials) &
       result(direct)
    complex(dp), intent(in) :: amplitudes(:)
    real(dp), intent(in) :: kernel(:)
    integer, intent(in) :: permutations(:, :), factorials(0:)

    complex(dp) :: direct
    integer, allocatable :: inverse(:), relative(:)
    integer :: column, position, rank, row

    allocate(inverse(size(permutations, 1)))
    allocate(relative(size(permutations, 1)))
    direct = cmplx(0.0_dp, 0.0_dp, dp)
    do row = 1, size(amplitudes)
      do position = 1, size(permutations, 1)
        inverse(permutations(position, row)) = position
      end do
      do column = 1, size(amplitudes)
        do position = 1, size(permutations, 1)
          relative(position) = inverse(permutations(position, column))
        end do
        rank = permutation_rank(relative, factorials)
        direct = direct+conjg(amplitudes(row))*kernel(rank)* &
             amplitudes(column)
      end do
    end do
  end function direct_colour_sum

  integer function permutation_rank(permutation, factorials) result(rank)
    integer, intent(in) :: permutation(:), factorials(0:)

    integer :: left, right, smaller

    rank = 1
    do left = 1, size(permutation)-1
      smaller = 0
      do right = left+1, size(permutation)
        if (permutation(right) < permutation(left)) smaller = smaller+1
      end do
      rank = rank+smaller*factorials(size(permutation)-left)
    end do
  end function permutation_rank

  integer function factorial(number) result(value)
    integer, intent(in) :: number

    integer :: factor

    value = 1
    do factor = 2, number
      if (value > huge(value)/factor) call fail('test factorial overflow')
      value = value*factor
    end do
  end function factorial

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) call fail(message)
  end subroutine require

end program test_trace_colour_matrix
