module adjoint_colour_matrix
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ampligluon_common, only: dp, fail, i64
  use symmetric_group_fft, only: symmetric_group_fft_t, &
       symmetric_group_fft_workspace_t
  use trace_colour_kernel, only: build_u3_trace_kernel
  implicit none
  private

  integer(i64), parameter :: number_of_colours = 3_i64
  real(dp), parameter :: kernel_imaginary_tolerance = &
       128.0_dp*epsilon(1.0_dp)
  real(dp), parameter :: kernel_symmetry_tolerance = 1.0e-10_dp

  ! For lexicographically ordered DDM permutations g,h in S_(N-2), the
  ! adjoint Gram matrix is the group convolution
  !
  !   C(g,h) = k(g^{-1} h).
  !
  ! The Fourier transform stores this dense matrix in (N-2)! coefficients.
  ! The nested commutators factor into left and right finite differences of
  ! the U(3) trace kernel.  Thus the DDM kernel needs only 2*(N-2) linear
  ! sweeps over S_(N-1), instead of (N-2)!*4^(N-2) trace pairings.
  type, public :: adjoint_colour_matrix_t
    private
    logical :: initialized = .false.
    logical :: use_colour_fft_value = .true.
    integer :: number_of_basis_elements_value = 0
    type(symmetric_group_fft_t), allocatable :: fourier_plan
    type(symmetric_group_fft_workspace_t), allocatable :: fourier_workspace
    real(dp), allocatable :: kernel_fourier(:)
    complex(dp), allocatable :: amplitude_fourier(:)
    real(dp), allocatable :: kernel_direct(:)
    integer, allocatable :: permutations_direct(:, :)
    integer(i64), allocatable :: factorials_direct(:)
    integer, allocatable :: inverse_direct(:)
    integer, allocatable :: relative_direct(:)
  contains
    procedure, public :: contract => contract_adjoint_colour_matrix
    procedure, public :: contract_reusing_workspace => &
         contract_adjoint_colour_matrix_reusing_workspace
    procedure, public :: initialize => initialize_adjoint_colour_matrix
    procedure, public :: number_of_basis_elements
    procedure, public :: number_of_stored_kernel_values
  end type adjoint_colour_matrix_t

contains

  subroutine initialize_adjoint_colour_matrix(this, total_gluons, permutations, &
       use_colour_fft)
    class(adjoint_colour_matrix_t), intent(inout) :: this
    integer, intent(in) :: total_gluons
    integer, allocatable, intent(out) :: permutations(:, :)
    logical, intent(in), optional :: use_colour_fft

    complex(dp), allocatable :: kernel_values(:)
    integer :: middle_length, order, position, trace_length
    integer(i64) :: number_of_basis_64
    integer(i64), allocatable :: factorials(:)
    integer, allocatable :: middle_permutations(:, :)
    integer(i64), allocatable :: adjoint_kernel(:)
    integer(i64) :: expected_diagonal
    logical :: enable_colour_fft

    if (total_gluons < 4) &
         call fail('gg scattering needs at least four total gluons')
    this%initialized = .false.
    this%use_colour_fft_value = .true.
    if (present(use_colour_fft)) this%use_colour_fft_value = use_colour_fft
    enable_colour_fft = this%use_colour_fft_value
    if (allocated(this%fourier_plan)) deallocate(this%fourier_plan)
    if (allocated(this%fourier_workspace)) deallocate(this%fourier_workspace)
    if (allocated(this%kernel_fourier)) deallocate(this%kernel_fourier)
    if (allocated(this%amplitude_fourier)) &
         deallocate(this%amplitude_fourier)
    if (allocated(this%kernel_direct)) deallocate(this%kernel_direct)
    if (allocated(this%permutations_direct)) &
         deallocate(this%permutations_direct)
    if (allocated(this%factorials_direct)) &
         deallocate(this%factorials_direct)
    if (allocated(this%inverse_direct)) deallocate(this%inverse_direct)
    if (allocated(this%relative_direct)) deallocate(this%relative_direct)

    middle_length = total_gluons-2
    trace_length = total_gluons-1
    if (trace_length >= bit_size(0_i64)-1) &
         call fail('too many gluons for permutation bit masks')
    if (middle_length >= bit_size(0)-1) &
         call fail('too many gluons for the commutator expansion')

    allocate(factorials(0:trace_length))
    factorials(0) = 1_i64
    do position = 1, trace_length
      if (factorials(position-1) > huge(0_i64)/int(position, i64)) &
           call fail('the number of permutations overflows a 64-bit integer')
      factorials(position) = factorials(position-1)*int(position, i64)
    end do
    number_of_basis_64 = factorials(middle_length)
    if (number_of_basis_64 > int(huge(0), i64) .or. &
         factorials(trace_length) > int(huge(0), i64)) &
         call fail('the number of permutations exceeds array-index capacity')
    this%number_of_basis_elements_value = int(number_of_basis_64)

    allocate(middle_permutations(middle_length, &
         this%number_of_basis_elements_value))
    middle_permutations(:, 1) = [(position, position=1, middle_length)]
    do order = 2, this%number_of_basis_elements_value
      middle_permutations(:, order) = middle_permutations(:, order-1)
      call next_permutation(middle_permutations(:, order))
    end do
    if (allocated(permutations)) deallocate(permutations)
    allocate(permutations(trace_length, this%number_of_basis_elements_value))
    permutations(1, :) = 1
    permutations(2:trace_length, :) = middle_permutations+1

    ! A commutator has no U(1) component, so the U(3) cycle kernel is exactly
    ! equivalent to the SU(3) Fierz kernel after both DDM chains are expanded.
    ! It avoids a separate 2^N singlet enumeration for every trace order.
    call build_adjoint_kernel(total_gluons, factorials, adjoint_kernel)
    expected_diagonal = (number_of_colours**2-1_i64)* &
         integer_power(2_i64*number_of_colours, middle_length)
    if (adjoint_kernel(1) /= expected_diagonal) &
         call fail('unexpected adjoint colour-matrix normalization')

    if (enable_colour_fft) then
      allocate(this%fourier_plan, this%fourier_workspace)
      call this%fourier_plan%initialize(middle_length)
      if (this%fourier_plan%group_order() /= &
           this%number_of_basis_elements_value) &
           call fail('adjoint kernel and Fourier group orders disagree')
      call this%fourier_plan%initialize_workspace(this%fourier_workspace)
      allocate(kernel_values(this%number_of_basis_elements_value))
      kernel_values = cmplx(real(adjoint_kernel, dp), 0.0_dp, dp)
      allocate(this%amplitude_fourier(this%number_of_basis_elements_value))
      call this%fourier_plan%forward_into(kernel_values, &
           this%amplitude_fourier, this%fourier_workspace)
      call store_real_kernel(this, this%amplitude_fourier)
    else
      allocate(this%kernel_direct(this%number_of_basis_elements_value))
      this%kernel_direct = real(adjoint_kernel, dp)
      call move_alloc(middle_permutations, this%permutations_direct)
      allocate(this%factorials_direct(0:middle_length))
      this%factorials_direct = factorials(0:middle_length)
      allocate(this%inverse_direct(middle_length), &
           this%relative_direct(middle_length))
    end if

    this%initialized = .true.
  end subroutine initialize_adjoint_colour_matrix

  subroutine build_adjoint_kernel(total_gluons, factorials, kernel)
    integer, intent(in) :: total_gluons
    integer(i64), intent(in) :: factorials(0:)
    integer(i64), allocatable, intent(out) :: kernel(:)

    integer :: allocation_status, generator, label, order, position, rank
    integer :: trace_length
    integer, allocatable :: permutation(:, :), relative(:)
    integer(i64), allocatable :: current(:), previous(:)
    integer(i64) :: magnitude_bound
    real(dp), allocatable :: trace_kernel(:)

    trace_length = total_gluons-1
    call build_u3_trace_kernel(total_gluons, trace_kernel)
    if (size(trace_kernel) /= int(factorials(trace_length))) &
         call fail('wrong trace kernel size in adjoint finite differences')
    allocate(permutation(trace_length, size(trace_kernel)), &
         stat=allocation_status)
    if (allocation_status /= 0) &
         call fail('not enough memory for adjoint kernel permutations')
    permutation(:, 1) = [(position, position=1, trace_length)]
    do order = 2, size(trace_kernel)
      permutation(:, order) = permutation(:, order-1)
      call next_permutation(permutation(:, order))
    end do
    allocate(current(size(trace_kernel)), previous(size(trace_kernel)), &
         stat=allocation_status)
    if (allocation_status /= 0) &
         call fail('not enough memory for adjoint finite differences')
    allocate(relative(trace_length), stat=allocation_status)
    if (allocation_status /= 0) &
         call fail('not enough memory for adjoint permutation workspace')
    ! build_u3_trace_kernel obtains each value from an exact checked i64 power.
    current = int(trace_kernel, i64)
    magnitude_bound = integer_power(number_of_colours, total_gluons)
    do generator = 1, 2*(total_gluons-2)
      if (magnitude_bound > shiftr(huge(magnitude_bound), 1)) &
           call fail('adjoint finite differences can overflow 64-bit integers')
      magnitude_bound = 2_i64*magnitude_bound
    end do

    ! Let s_i=[i,1,...,i-1,i+1,...,L] in S_L.  The canonical DDM trace
    ! expansion is F=(1-s_2)(1-s_3)...(1-s_L).  Applying F on both sides
    ! of the trace Gram kernel gives these exact finite differences:
    !
    !   left:  f(a) <- f(a)-f(s_i^{-1} a)
    !   right: f(a) <- f(a)-f(a s_i).
    !
    ! Each sweep reads the previous complete group function; updating in
    ! place would be incorrect because s_i is a cycle rather than an involution.
    do generator = trace_length, 2, -1
      previous = current
      do order = 1, size(current)
        do position = 1, trace_length
          label = permutation(position, order)
          if (label < generator) then
            relative(position) = label+1
          else if (label == generator) then
            relative(position) = 1
          else
            relative(position) = label
          end if
        end do
        rank = permutation_rank(relative, factorials)
        current(order) = previous(order)-previous(rank)
      end do
    end do
    do generator = trace_length, 2, -1
      previous = current
      do order = 1, size(current)
        relative(1) = permutation(generator, order)
        relative(2:generator) = permutation(1:generator-1, order)
        if (generator < trace_length) &
             relative(generator+1:) = permutation(generator+1:, order)
        rank = permutation_rank(relative, factorials)
        current(order) = previous(order)-previous(rank)
      end do
    end do

    ! The lexicographic prefix of S_L consists precisely of [1,sigma+1]
    ! for sigma in S_(L-1), which is the public DDM basis order.
    allocate(kernel(int(factorials(trace_length-1))), stat=allocation_status)
    if (allocation_status /= 0) &
         call fail('not enough memory for adjoint colour kernel')
    kernel = current(1:size(kernel))
  end subroutine build_adjoint_kernel

  subroutine contract_adjoint_colour_matrix(this, amplitudes, result)
    class(adjoint_colour_matrix_t), intent(in) :: this
    complex(dp), intent(in) :: amplitudes(:)
    real(dp), intent(out) :: result

    complex(dp), allocatable :: transformed_amplitudes(:)

    if (.not. this%initialized) &
         call fail('adjoint colour matrix is not initialized')
    if (size(amplitudes) /= this%number_of_basis_elements_value) &
         call fail('wrong number of adjoint-basis amplitudes')

    if (this%use_colour_fft_value) then
      if (.not. allocated(this%fourier_plan)) &
           call fail('adjoint colour Fourier plan is not initialized')
      call this%fourier_plan%forward(amplitudes, transformed_amplitudes)
      call contract_fourier_blocks(this, transformed_amplitudes, result)
    else
      call contract_direct(this, amplitudes, result)
    end if
  end subroutine contract_adjoint_colour_matrix

  subroutine contract_adjoint_colour_matrix_reusing_workspace( &
       this, amplitudes, result)
    class(adjoint_colour_matrix_t), intent(inout) :: this
    complex(dp), intent(in) :: amplitudes(:)
    real(dp), intent(out) :: result

    if (.not. this%initialized) &
         call fail('adjoint colour matrix is not initialized')
    if (size(amplitudes) /= this%number_of_basis_elements_value) &
         call fail('wrong number of adjoint-basis amplitudes')

    if (this%use_colour_fft_value) then
      if (.not. allocated(this%fourier_plan) .or. &
           .not. allocated(this%fourier_workspace)) &
           call fail('adjoint colour Fourier workspace is not initialized')
      call this%fourier_plan%forward_into(amplitudes, &
           this%amplitude_fourier, this%fourier_workspace)
      call contract_fourier_blocks(this, this%amplitude_fourier, result)
    else
      if (.not. allocated(this%inverse_direct) .or. &
           .not. allocated(this%relative_direct)) &
           call fail('direct adjoint colour workspace is not initialized')
      call contract_direct_with_workspace(this, amplitudes, result, &
           this%inverse_direct, this%relative_direct)
    end if
  end subroutine contract_adjoint_colour_matrix_reusing_workspace

  subroutine contract_direct(this, amplitudes, result)
    class(adjoint_colour_matrix_t), intent(in) :: this
    complex(dp), intent(in) :: amplitudes(:)
    real(dp), intent(out) :: result

    integer, allocatable :: inverse(:), relative(:)

    if (.not. allocated(this%kernel_direct) .or. &
         .not. allocated(this%permutations_direct) .or. &
         .not. allocated(this%factorials_direct)) &
         call fail('direct adjoint colour matrix is not initialized')
    allocate(inverse(size(this%permutations_direct, 1)), &
         relative(size(this%permutations_direct, 1)))
    call contract_direct_with_workspace(this, amplitudes, result, inverse, relative)
  end subroutine contract_direct

  subroutine contract_direct_with_workspace(this, amplitudes, result, inverse, &
       relative)
    class(adjoint_colour_matrix_t), intent(in) :: this
    complex(dp), intent(in) :: amplitudes(:)
    real(dp), intent(out) :: result
    integer, intent(inout) :: inverse(:), relative(:)

    integer :: column, position, rank, row
    real(dp) :: total

    if (size(inverse) /= size(this%permutations_direct, 1) .or. &
         size(relative) /= size(this%permutations_direct, 1)) &
         call fail('direct adjoint colour workspace has the wrong size')

    ! Every diagonal element is k(e).  Since the Gram matrix is real
    ! symmetric, each off-diagonal amplitude overlap is needed only once.
    total = 0.0_dp
    do row = 1, this%number_of_basis_elements_value
      total = total+this%kernel_direct(1)*( &
           real(amplitudes(row), dp)**2+aimag(amplitudes(row))**2)
    end do
    do row = 1, this%number_of_basis_elements_value-1
      do position = 1, size(this%permutations_direct, 1)
        inverse(this%permutations_direct(position, row)) = position
      end do
      do column = row+1, this%number_of_basis_elements_value
        do position = 1, size(this%permutations_direct, 1)
          relative(position) = inverse( &
               this%permutations_direct(position, column))
        end do
        rank = permutation_rank(relative, this%factorials_direct)
        total = total+2.0_dp*this%kernel_direct(rank)*real( &
             conjg(amplitudes(row))*amplitudes(column), dp)
      end do
    end do
    if (.not. ieee_is_finite(total)) &
         call fail('direct adjoint colour contraction produced a non-finite value')
    result = total
  end subroutine contract_direct_with_workspace

  subroutine contract_fourier_blocks(this, transformed_amplitudes, result)
    class(adjoint_colour_matrix_t), intent(in) :: this
    complex(dp), intent(in) :: transformed_amplitudes(:)
    real(dp), intent(out) :: result

    complex(dp) :: first_amplitude1, first_amplitude2
    complex(dp) :: first_amplitude3, first_amplitude4, other_amplitude
    real(dp) :: block_sum, column_norm, overlaps1(4), overlaps2(4)
    real(dp) :: overlaps3(4), overlaps4(4), total
    integer :: base(4), block, column, column_base, dimension, first
    integer :: lane, offset, other, row

    if (size(transformed_amplitudes) /= &
         this%number_of_basis_elements_value) &
         call fail('wrong transformed adjoint amplitude size')

    ! With F_lambda = sum_g f(g) rho_lambda(g), right convolution by
    ! k(g^{-1} h) transforms as A_lambda K_lambda.  Plancherel gives
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
               real(transformed_amplitudes(column_base+row), dp)**2+ &
               aimag(transformed_amplitudes(column_base+row))**2
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
                 transformed_amplitudes(base(1)+row))
            first_amplitude2 = conjg( &
                 transformed_amplitudes(base(2)+row))
            first_amplitude3 = conjg( &
                 transformed_amplitudes(base(3)+row))
            first_amplitude4 = conjg( &
                 transformed_amplitudes(base(4)+row))
            do lane = 0, 3
              other_amplitude = transformed_amplitudes( &
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
            other_amplitude = transformed_amplitudes( &
                 offset+(other-1)*dimension+row)
            do first = 1, 4
              overlaps1(first) = overlaps1(first)+real(conjg( &
                   transformed_amplitudes(base(first)+row))* &
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
                   transformed_amplitudes(base(first)+row))* &
                   transformed_amplitudes(base(other)+row), dp)
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
                 transformed_amplitudes(column_base+row))* &
                 transformed_amplitudes( &
                 offset+(other-1)*dimension+row), dp)
          end do
          block_sum = block_sum+2.0_dp*this%kernel_fourier( &
               offset+(other-1)*dimension+first-1)*overlaps1(1)
        end do
      end do
      total = total+real(dimension, dp)*block_sum
    end do
    if (.not. ieee_is_finite(total)) &
         call fail('adjoint colour contraction produced a non-finite value')
    result = total/real(this%number_of_basis_elements_value, dp)
  end subroutine contract_fourier_blocks

  subroutine store_real_kernel(this, transformed_kernel)
    class(adjoint_colour_matrix_t), intent(inout) :: this
    complex(dp), intent(in) :: transformed_kernel(:)

    real(dp) :: block_scale, first_value, second_value
    integer :: block, column, dimension, first_index, offset, row
    integer :: second_index

    if (size(transformed_kernel) /= this%number_of_basis_elements_value) &
         call fail('transformed adjoint colour kernel has the wrong size')
    if (.not. all(ieee_is_finite(real(transformed_kernel, dp))) .or. &
         .not. all(ieee_is_finite(aimag(transformed_kernel)))) &
         call fail('transformed adjoint colour kernel is not finite')
    block_scale = max(1.0_dp, maxval(abs(transformed_kernel)))
    if (maxval(abs(aimag(transformed_kernel))) > &
         kernel_imaginary_tolerance*block_scale) &
         call fail('transformed real adjoint kernel has an imaginary part')

    allocate(this%kernel_fourier(this%number_of_basis_elements_value))
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
               kernel_symmetry_tolerance*block_scale) &
               call fail('transformed adjoint colour kernel is not symmetric')
          this%kernel_fourier(first_index) = &
               0.5_dp*(first_value+second_value)
        end do
      end do
    end do
  end subroutine store_real_kernel

  integer function number_of_basis_elements(this) result(number)
    class(adjoint_colour_matrix_t), intent(in) :: this

    if (.not. this%initialized) &
         call fail('adjoint colour matrix is not initialized')
    number = this%number_of_basis_elements_value
  end function number_of_basis_elements

  integer function number_of_stored_kernel_values(this) result(number)
    class(adjoint_colour_matrix_t), intent(in) :: this

    if (.not. this%initialized) &
         call fail('adjoint colour matrix is not initialized')
    if (this%use_colour_fft_value) then
      if (.not. allocated(this%kernel_fourier)) &
           call fail('adjoint Fourier colour kernel is not initialized')
      number = size(this%kernel_fourier)
    else
      if (.not. allocated(this%kernel_direct)) &
           call fail('direct adjoint colour kernel is not initialized')
      number = size(this%kernel_direct)
    end if
  end function number_of_stored_kernel_values

  integer function permutation_rank(permutation, factorials) result(rank)
    integer, intent(in) :: permutation(:)
    integer(i64), intent(in) :: factorials(0:)

    integer :: label, position, smaller
    integer(i64) :: lower_mask, rank_64, used_mask

    rank_64 = 1_i64
    used_mask = 0_i64
    do position = 1, size(permutation)-1
      label = permutation(position)
      lower_mask = shiftl(1_i64, label-1)-1_i64
      smaller = label-1-popcnt(iand(used_mask, lower_mask))
      rank_64 = rank_64+int(smaller, i64)* &
           factorials(size(permutation)-position)
      used_mask = ibset(used_mask, label-1)
    end do
    if (rank_64 > int(huge(0), i64)) &
         call fail('permutation rank exceeds array-index capacity')
    rank = int(rank_64)
  end function permutation_rank

  subroutine next_permutation(permutation)
    integer, intent(inout) :: permutation(:)

    integer :: left, right, temporary

    left = size(permutation)-1
    do while (left >= 1)
      if (permutation(left) < permutation(left+1)) exit
      left = left-1
    end do
    if (left < 1) then
      permutation = permutation(size(permutation):1:-1)
      return
    end if
    right = size(permutation)
    do while (permutation(right) < permutation(left))
      right = right-1
    end do
    temporary = permutation(left)
    permutation(left) = permutation(right)
    permutation(right) = temporary
    permutation(left+1:) = permutation(size(permutation):left+1:-1)
  end subroutine next_permutation

  integer(i64) function integer_power(base, exponent) result(value)
    integer(i64), intent(in) :: base
    integer, intent(in) :: exponent

    integer :: power

    if (exponent < 0) call fail('negative exponent in exact colour factor')
    value = 1_i64
    do power = 1, exponent
      if (value > huge(value)/base) &
           call fail('colour factor power overflows a 64-bit integer')
      value = value*base
    end do
  end function integer_power

end module adjoint_colour_matrix
