module trace_colour_matrix
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ampligluon_common, only: dp, fail, i64
  use symmetric_group_fft, only: symmetric_group_fft_t, &
       symmetric_group_fft_workspace_t
  use trace_colour_kernel, only: build_u3_trace_kernel
  implicit none
  private

  ! For lexicographically ordered permutations g,h in S_(N-1), the trace
  ! Gram matrix is the group convolution
  !
  !   C(g,h) = k(g^{-1} h).
  !
  ! The Fourier transform turns this dense matrix into one small matrix for
  ! every Young irrep.  Both the transformed kernel and an amplitude transform
  ! contain sum_lambda d_lambda**2 = (N-1)! coefficients.
  real(dp), parameter :: kernel_imaginary_tolerance = &
       128.0_dp*epsilon(1.0_dp)
  real(dp), parameter :: kernel_symmetry_tolerance_floor = 1.0e-10_dp
  real(dp), parameter :: kernel_symmetry_roundoff_factor = 512.0_dp

  type, public :: trace_colour_matrix_t
    private
    logical :: initialized = .false.
    logical :: use_colour_fft_value = .true.
    integer :: number_of_orders_value = 0
    type(symmetric_group_fft_t), allocatable :: fourier_plan
    type(symmetric_group_fft_workspace_t), allocatable :: fourier_workspace
    real(dp), allocatable :: kernel_fourier(:)
    real(dp), allocatable :: kernel_direct(:)
    complex(dp), allocatable :: amplitude_fourier(:)
    integer(i64), allocatable :: factorials(:)
    integer, allocatable :: permutations(:, :)
  contains
    procedure, public :: contract => contract_trace_colour_matrix
    procedure, public :: initialize => initialize_trace_colour_matrix
    procedure, public :: number_of_orders
    procedure, public :: number_of_stored_kernel_values
  end type trace_colour_matrix_t

contains

  subroutine initialize_trace_colour_matrix(this, total_gluons, permutations, &
       use_colour_fft)
    class(trace_colour_matrix_t), intent(inout) :: this
    integer, intent(in) :: total_gluons
    integer, allocatable, intent(out), optional :: permutations(:, :)
    logical, intent(in), optional :: use_colour_fft

    complex(dp), allocatable :: kernel_values(:)
    real(dp), allocatable :: kernel(:)
    integer(i64) :: number_of_orders_64
    integer :: length, position

    if (total_gluons < 4) &
         call fail('gg scattering needs at least four total gluons')
    this%initialized = .false.
    this%use_colour_fft_value = .true.
    if (present(use_colour_fft)) &
         this%use_colour_fft_value = use_colour_fft
    if (allocated(this%fourier_plan)) deallocate(this%fourier_plan)
    if (allocated(this%fourier_workspace)) &
         deallocate(this%fourier_workspace)
    if (allocated(this%kernel_fourier)) deallocate(this%kernel_fourier)
    if (allocated(this%kernel_direct)) deallocate(this%kernel_direct)
    if (allocated(this%amplitude_fourier)) &
         deallocate(this%amplitude_fourier)
    if (allocated(this%factorials)) deallocate(this%factorials)
    if (allocated(this%permutations)) deallocate(this%permutations)
    length = total_gluons-1
    number_of_orders_64 = factorial_i64(length)
    if (number_of_orders_64 > int(huge(0), i64)) &
         call fail('the number of colour orders exceeds array-index capacity')
    this%number_of_orders_value = int(number_of_orders_64)

    call build_u3_trace_kernel(total_gluons, kernel)
    if (size(kernel) /= this%number_of_orders_value) &
         call fail('trace colour kernel has the wrong group order')
    if (.not. all(ieee_is_finite(kernel))) &
         call fail('trace colour kernel is not finite')

    if (this%use_colour_fft_value) then
      if (present(permutations)) call build_permutation_table(length, &
           this%number_of_orders_value, permutations)
      allocate(this%fourier_plan, this%fourier_workspace)
      call this%fourier_plan%initialize(length)
      if (this%fourier_plan%group_order() /= this%number_of_orders_value) &
           call fail('colour kernel and Fourier group orders disagree')
      call this%fourier_plan%initialize_workspace(this%fourier_workspace)
      allocate(kernel_values(this%number_of_orders_value))
      kernel_values = cmplx(kernel, 0.0_dp, dp)
      deallocate(kernel)
      allocate(this%amplitude_fourier(this%number_of_orders_value))
      call this%fourier_plan%forward_into(kernel_values, &
           this%amplitude_fourier, this%fourier_workspace)
      call store_real_kernel(this, this%amplitude_fourier)
    else
      allocate(this%factorials(0:length))
      this%factorials(0) = 1_i64
      do position = 1, length
        this%factorials(position) = this%factorials(position-1)* &
             int(position, i64)
      end do
      if (this%factorials(length) /= number_of_orders_64) &
           call fail('direct colour factorials have the wrong group order')
      call build_permutation_table(length, this%number_of_orders_value, &
           this%permutations)
      if (present(permutations)) permutations = this%permutations
      call move_alloc(kernel, this%kernel_direct)
    end if

    this%initialized = .true.
  end subroutine initialize_trace_colour_matrix

  subroutine contract_trace_colour_matrix(this, amplitudes, result)
    class(trace_colour_matrix_t), intent(inout) :: this
    complex(dp), intent(in) :: amplitudes(:)
    real(dp), intent(out) :: result

    complex(dp) :: first_amplitude1, first_amplitude2
    complex(dp) :: first_amplitude3, first_amplitude4, other_amplitude
    real(dp) :: block_sum, column_norm, overlaps1(4), overlaps2(4)
    real(dp) :: overlaps3(4), overlaps4(4), total
    integer :: base(4), block, column, column_base, dimension, first
    integer :: lane, offset, other, row

    if (.not. this%initialized) &
         call fail('colour contraction is not initialized')
    if (size(amplitudes) /= this%number_of_orders_value) &
         call fail('wrong number of colour-ordered amplitudes')

    if (.not. this%use_colour_fft_value) then
      call contract_trace_colour_matrix_direct(this, amplitudes, result)
      return
    end if
    if (.not. allocated(this%fourier_plan) .or. &
         .not. allocated(this%fourier_workspace) .or. &
         .not. allocated(this%amplitude_fourier)) &
         call fail('trace Fourier colour contraction has no workspace')

    call this%fourier_plan%forward_into(amplitudes, &
         this%amplitude_fourier, this%fourier_workspace)

    ! With F_lambda = sum_g f(g) rho_lambda(g), right convolution by
    ! k(g^{-1} h) transforms as A_lambda K_lambda.  Plancherel therefore
    ! gives
    !
    !   A^dagger C A = (1/|G|) sum_lambda d_lambda
    !                       Tr(A_lambda^dagger A_lambda K_lambda).
    total = 0.0_dp
    do block = 1, this%fourier_plan%number_of_blocks()
      dimension = this%fourier_plan%block_dimension(block)
      offset = this%fourier_plan%block_offset(block)
      block_sum = 0.0_dp

      ! Accumulate the diagonal separately.  The off-diagonal microkernel
      ! below then reuses every streamed column for two column overlaps.
      do column = 1, dimension
        column_base = offset+(column-1)*dimension
        column_norm = 0.0_dp
        do row = 0, dimension-1
          column_norm = column_norm+ &
               real(this%amplitude_fourier(column_base+row), dp)**2+ &
               aimag(this%amplitude_fourier(column_base+row))**2
        end do
        block_sum = block_sum+this%kernel_fourier( &
             column_base+column-1)*column_norm
      end do

      ! K_lambda is real symmetric.  Pair the (column,other) and
      ! (other,column) terms so each complex column overlap is evaluated
      ! once and the quadratic form stays explicitly real:
      !
      !   Tr(A^dagger A K) = sum_i K_ii ||A_i||^2
      !       + 2 sum_(i<j) K_ij Re(A_i^dagger A_j).
      !
      ! Work on four first columns and four other columns at a time.  The
      ! sixteen independent reductions expose instruction-level parallelism
      ! and reuse each loaded amplitude in four column overlaps.
      column = 1
      do while (column+3 <= dimension)
        do first = 1, 4
          base(first) = offset+(column+first-2)*dimension
        end do
        other = column+4
        do while (other+3 <= dimension)
          overlaps1 = 0.0_dp
          overlaps2 = 0.0_dp
          overlaps3 = 0.0_dp
          overlaps4 = 0.0_dp
          do row = 0, dimension-1
            first_amplitude1 = conjg( &
                 this%amplitude_fourier(base(1)+row))
            first_amplitude2 = conjg( &
                 this%amplitude_fourier(base(2)+row))
            first_amplitude3 = conjg( &
                 this%amplitude_fourier(base(3)+row))
            first_amplitude4 = conjg( &
                 this%amplitude_fourier(base(4)+row))
            do lane = 0, 3
              other_amplitude = this%amplitude_fourier( &
                   offset+(other+lane-1)*dimension+row)
              overlaps1(lane+1) = overlaps1(lane+1)+ &
                   real(first_amplitude1*other_amplitude, dp)
              overlaps2(lane+1) = overlaps2(lane+1)+ &
                   real(first_amplitude2*other_amplitude, dp)
              overlaps3(lane+1) = overlaps3(lane+1)+ &
                   real(first_amplitude3*other_amplitude, dp)
              overlaps4(lane+1) = overlaps4(lane+1)+ &
                   real(first_amplitude4*other_amplitude, dp)
            end do
          end do
          do lane = 0, 3
            block_sum = block_sum+2.0_dp*this%kernel_fourier( &
                 offset+(other+lane-1)*dimension+column-1)* &
                 overlaps1(lane+1)
            block_sum = block_sum+2.0_dp*this%kernel_fourier( &
                 offset+(other+lane-1)*dimension+column)* &
                 overlaps2(lane+1)
            block_sum = block_sum+2.0_dp*this%kernel_fourier( &
                 offset+(other+lane-1)*dimension+column+1)* &
                 overlaps3(lane+1)
            block_sum = block_sum+2.0_dp*this%kernel_fourier( &
                 offset+(other+lane-1)*dimension+column+2)* &
                 overlaps4(lane+1)
          end do
          other = other+4
        end do
        do while (other <= dimension)
          overlaps1 = 0.0_dp
          do row = 0, dimension-1
            other_amplitude = this%amplitude_fourier( &
                 offset+(other-1)*dimension+row)
            do first = 1, 4
              overlaps1(first) = overlaps1(first)+real(conjg( &
                   this%amplitude_fourier(base(first)+row))* &
                   other_amplitude, dp)
            end do
          end do
          do first = 1, 4
            block_sum = block_sum+2.0_dp*this%kernel_fourier( &
                 offset+(other-1)*dimension+column+first-2)* &
                 overlaps1(first)
          end do
          other = other+1
        end do

        ! Finish the six overlaps within the four first columns.
        do first = 1, 3
          do other = first+1, 4
            overlaps1(1) = 0.0_dp
            do row = 0, dimension-1
              overlaps1(1) = overlaps1(1)+real(conjg( &
                   this%amplitude_fourier(base(first)+row))* &
                   this%amplitude_fourier(base(other)+row), dp)
            end do
            block_sum = block_sum+2.0_dp*this%kernel_fourier( &
                 base(other)+column+first-2)*overlaps1(1)
          end do
        end do
        column = column+4
      end do

      ! At most three first columns remain.  Their overlaps with earlier
      ! tiles were already accumulated, so only their upper triangle remains.
      do first = column, dimension-1
        column_base = offset+(first-1)*dimension
        do other = first+1, dimension
          overlaps1(1) = 0.0_dp
          do row = 0, dimension-1
            overlaps1(1) = overlaps1(1)+real(conjg( &
                 this%amplitude_fourier(column_base+row))* &
                 this%amplitude_fourier( &
                 offset+(other-1)*dimension+row), dp)
          end do
          block_sum = block_sum+2.0_dp*this%kernel_fourier( &
               offset+(other-1)*dimension+first-1)*overlaps1(1)
        end do
      end do
      total = total+real(dimension, dp)*block_sum
    end do
    if (.not. ieee_is_finite(total)) &
         call fail('trace colour contraction produced a non-finite value')
    result = total/real(this%number_of_orders_value, dp)
  end subroutine contract_trace_colour_matrix

  subroutine contract_trace_colour_matrix_direct(this, amplitudes, result)
    class(trace_colour_matrix_t), intent(in) :: this
    complex(dp), intent(in) :: amplitudes(:)
    real(dp), intent(out) :: result

    integer, allocatable :: inverse(:), relative(:)
    integer :: column, position, rank, row
    real(dp) :: column_imaginary, column_real, row_imaginary, row_real, total

    if (.not. allocated(this%kernel_direct) .or. &
         .not. allocated(this%factorials) .or. &
         .not. allocated(this%permutations)) &
         call fail('direct trace colour contraction has no kernel data')
    if (size(this%kernel_direct) /= this%number_of_orders_value .or. &
         size(this%permutations, 2) /= this%number_of_orders_value) &
         call fail('direct trace colour contraction has inconsistent storage')
    allocate(inverse(size(this%permutations, 1)))
    allocate(relative(size(this%permutations, 1)))

    ! The raw kernel is the identity row of C(g,h)=k(g^{-1}h).  Construct
    ! each other entry by mapping the lexicographic column permutation through
    ! the inverse row permutation.  Inversion symmetry makes C real symmetric,
    ! so only the strict upper triangle is needed.
    total = 0.0_dp
    do row = 1, this%number_of_orders_value
      row_real = real(amplitudes(row), dp)
      row_imaginary = aimag(amplitudes(row))
      total = total+this%kernel_direct(1)* &
           (row_real*row_real+row_imaginary*row_imaginary)
      do position = 1, size(this%permutations, 1)
        inverse(this%permutations(position, row)) = position
      end do
      do column = row+1, this%number_of_orders_value
        do position = 1, size(this%permutations, 1)
          relative(position) = inverse(this%permutations(position, column))
        end do
        rank = permutation_rank(relative, this%factorials)
        column_real = real(amplitudes(column), dp)
        column_imaginary = aimag(amplitudes(column))
        total = total+2.0_dp*this%kernel_direct(rank)* &
             (row_real*column_real+row_imaginary*column_imaginary)
      end do
    end do
    if (.not. ieee_is_finite(total)) &
         call fail('direct trace colour contraction produced a non-finite value')
    result = total
  end subroutine contract_trace_colour_matrix_direct

  subroutine store_real_kernel(this, transformed_kernel)
    class(trace_colour_matrix_t), intent(inout) :: this
    complex(dp), intent(in) :: transformed_kernel(:)

    real(dp) :: block_scale, first_value, second_value, symmetry_tolerance
    integer :: block, column, dimension, first_index, offset, row
    integer :: second_index

    if (size(transformed_kernel) /= this%number_of_orders_value) &
         call fail('transformed trace colour kernel has the wrong size')
    if (.not. all(ieee_is_finite(real(transformed_kernel, dp))) .or. &
         .not. all(ieee_is_finite(aimag(transformed_kernel)))) &
         call fail('transformed trace colour kernel is not finite')
    block_scale = max(1.0_dp, maxval(abs(transformed_kernel)))
    if (maxval(abs(aimag(transformed_kernel))) > &
         kernel_imaginary_tolerance*block_scale) &
         call fail('transformed real trace kernel has an imaginary part')

    ! Exact inversion symmetry makes every real Fourier block symmetric.
    ! Cancellation-dominated blocks can nevertheless retain an antisymmetric
    ! binary64 residue whose RMS accumulation envelope grows as sqrt(|G|).
    ! Keep the established small-group floor; 512 is the smallest power-of-two
    ! guard factor that covers the multi-stage transform through supported S_11.
    symmetry_tolerance = max(kernel_symmetry_tolerance_floor, &
         kernel_symmetry_roundoff_factor*epsilon(1.0_dp)* &
         sqrt(real(this%number_of_orders_value, dp)))
    allocate(this%kernel_fourier(this%number_of_orders_value))
    do block = 1, this%fourier_plan%number_of_blocks()
      dimension = this%fourier_plan%block_dimension(block)
      offset = this%fourier_plan%block_offset(block)
      block_scale = max(1.0_dp, maxval(abs(real(transformed_kernel( &
           offset:offset+dimension**2-1), dp))))
      do column = 1, dimension
        do row = 1, dimension
          first_index = offset+(column-1)*dimension+row-1
          second_index = offset+(row-1)*dimension+column-1
          first_value = real(transformed_kernel(first_index), dp)
          second_value = real(transformed_kernel(second_index), dp)
          if (abs(first_value-second_value) > &
               symmetry_tolerance*block_scale) &
               call fail('transformed trace colour kernel is not symmetric')
          ! Averaging makes the validated symmetry exact in stored arithmetic.
          this%kernel_fourier(first_index) = &
               0.5_dp*(first_value+second_value)
        end do
      end do
    end do
  end subroutine store_real_kernel

  integer function number_of_orders(this) result(number)
    class(trace_colour_matrix_t), intent(in) :: this

    if (.not. this%initialized) &
         call fail('colour contraction is not initialized')
    number = this%number_of_orders_value
  end function number_of_orders

  integer function number_of_stored_kernel_values(this) result(number)
    class(trace_colour_matrix_t), intent(in) :: this

    if (.not. this%initialized) &
         call fail('colour contraction is not initialized')
    if (this%use_colour_fft_value) then
      number = size(this%kernel_fourier)
    else
      number = size(this%kernel_direct)
    end if
  end function number_of_stored_kernel_values

  subroutine build_permutation_table(length, number_of_orders, permutations)
    integer, intent(in) :: length, number_of_orders
    integer, allocatable, intent(out) :: permutations(:, :)

    integer :: order, position

    allocate(permutations(length, number_of_orders))
    permutations(:, 1) = [(position, position=1, length)]
    do order = 2, number_of_orders
      permutations(:, order) = permutations(:, order-1)
      call next_permutation(permutations(:, order))
    end do
  end subroutine build_permutation_table

  integer function permutation_rank(permutation, factorials) result(rank)
    integer, intent(in) :: permutation(:)
    integer(i64), intent(in) :: factorials(0:)

    integer :: left, right, smaller
    integer(i64) :: rank_64

    if (ubound(factorials, 1) < size(permutation)) &
         call fail('direct trace colour rank has too few factorials')
    rank_64 = 1_i64
    do left = 1, size(permutation)-1
      smaller = 0
      do right = left+1, size(permutation)
        if (permutation(right) < permutation(left)) smaller = smaller+1
      end do
      rank_64 = rank_64+int(smaller, i64)* &
           factorials(size(permutation)-left)
    end do
    if (rank_64 > int(huge(rank), i64)) &
         call fail('direct trace colour rank exceeds array-index capacity')
    rank = int(rank_64)
  end function permutation_rank

  integer(i64) function factorial_i64(number) result(value)
    integer, intent(in) :: number

    integer :: factor

    if (number < 0) call fail('negative colour-order factorial')
    value = 1_i64
    do factor = 2, number
      if (value > huge(value)/int(factor, i64)) &
           call fail('the number of colour orders overflows a 64-bit integer')
      value = value*int(factor, i64)
    end do
  end function factorial_i64

  subroutine next_permutation(permutation)
    integer, intent(inout) :: permutation(:)

    integer :: left, right, temporary

    left = size(permutation)-1
    do while (left >= 1)
      if (permutation(left) < permutation(left+1)) exit
      left = left-1
    end do
    if (left < 1) &
         call fail('colour-order permutation sequence wrapped unexpectedly')
    right = size(permutation)
    do while (permutation(right) < permutation(left))
      right = right-1
    end do
    temporary = permutation(left)
    permutation(left) = permutation(right)
    permutation(right) = temporary
    permutation(left+1:) = permutation(size(permutation):left+1:-1)
  end subroutine next_permutation

end module trace_colour_matrix
