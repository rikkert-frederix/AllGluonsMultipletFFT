program test_symmetric_group_fft
  use ampligluon_common, only: dp, fail
  use symmetric_group_fft, only: symmetric_group_fft_t, &
       symmetric_group_fft_workspace_t
  implicit none

  integer :: degree

  call check_known_layouts()
  do degree = 1, 6
    call check_identity_impulse(degree)
    call check_direct_transform(degree)
    call check_representation_relations(degree)
    call check_reusable_workspace(degree)
  end do
  call check_parseval_transform(7)
  call check_parseval_transform(8)
  call check_parseval_transform(9)
  call check_selected_impulse(7, 1379)
  call check_selected_impulse(8, 17011)
  call check_selected_impulse(9, 271829)
  write(*, '(a)') 'Symmetric-group FFT regression: PASS'

contains

  subroutine check_known_layouts()
    call check_layout(1, [1])
    call check_layout(2, [1, 1])
    call check_layout(3, [1, 2, 1])
    call check_layout(4, [1, 3, 2, 3, 1])
    call check_layout(5, [1, 4, 5, 6, 5, 4, 1])
    call check_layout(6, [1, 5, 9, 10, 5, 16, 10, 5, 9, 5, 1])
  end subroutine check_known_layouts

  subroutine check_layout(degree, expected_dimensions)
    integer, intent(in) :: degree, expected_dimensions(:)

    type(symmetric_group_fft_t) :: plan
    integer :: block, coefficient_count, shape_sum
    integer, allocatable :: dimensions(:), offsets(:), shape(:)

    call plan%initialize(degree)
    call require(plan%degree() == degree, 'wrong FFT degree')
    call require(plan%group_order() == factorial(degree), &
         'wrong symmetric-group order')
    call require(plan%number_of_blocks() == size(expected_dimensions), &
         'wrong number of partitions')
    call plan%block_layout(offsets, dimensions)
    call require(all(dimensions == expected_dimensions), &
         'wrong Young-irrep dimensions')

    coefficient_count = 0
    do block = 1, plan%number_of_blocks()
      call require(offsets(block) == coefficient_count+1, &
           'Fourier blocks are not contiguous')
      call require(plan%block_offset(block) == offsets(block), &
           'scalar and vector block offsets disagree')
      call require(plan%block_dimension(block) == dimensions(block), &
           'scalar and vector block dimensions disagree')
      coefficient_count = coefficient_count+dimensions(block)**2
      shape = plan%block_shape(block)
      shape_sum = sum(shape)
      call require(shape_sum == degree, 'partition has the wrong degree')
      call require(all(shape >= 0), 'partition has a negative row length')
      if (degree > 1) call require(all(shape(1:degree-1) >= &
           shape(2:degree)), 'partition rows are not ordered')
    end do
    call require(coefficient_count == factorial(degree), &
         'Fourier block sizes do not sum to the group order')
  end subroutine check_layout

  subroutine check_identity_impulse(degree)
    integer, intent(in) :: degree

    type(symmetric_group_fft_t) :: plan
    complex(dp), allocatable :: coefficients(:), values(:)
    integer :: block, column, dimension, offset, row

    call plan%initialize(degree)
    allocate(values(plan%group_order()))
    values = cmplx(0.0_dp, 0.0_dp, dp)
    values(1) = cmplx(1.0_dp, 0.0_dp, dp)
    call plan%forward(values, coefficients)
    do block = 1, plan%number_of_blocks()
      dimension = plan%block_dimension(block)
      offset = plan%block_offset(block)
      do column = 1, dimension
        do row = 1, dimension
          if (row == column) then
            call require_close(coefficients(offset+(column-1)*dimension+ &
                 row-1), cmplx(1.0_dp, 0.0_dp, dp), &
                 'identity impulse has a non-unit diagonal')
          else
            call require_close(coefficients(offset+(column-1)*dimension+ &
                 row-1), cmplx(0.0_dp, 0.0_dp, dp), &
                 'identity impulse has a nonzero off-diagonal entry')
          end if
        end do
      end do
    end do
  end subroutine check_identity_impulse

  subroutine check_direct_transform(degree)
    integer, intent(in) :: degree

    type(symmetric_group_fft_t) :: plan
    complex(dp), allocatable :: direct(:), fast(:), values(:)
    complex(dp) :: sign_transform
    real(dp), allocatable :: representation(:, :)
    real(dp) :: fourier_norm, input_norm, scale
    integer :: block, dimension, offset, order, permutation_index
    integer, allocatable :: permutation(:)

    call plan%initialize(degree)
    order = plan%group_order()
    allocate(values(order), direct(order), permutation(degree))
    do permutation_index = 1, order
      values(permutation_index) = cmplx( &
           sin(0.37_dp*real(permutation_index, dp))+ &
           0.013_dp*real(permutation_index, dp), &
           cos(0.23_dp*real(permutation_index, dp))- &
           0.007_dp*real(permutation_index, dp), dp)
    end do

    call plan%forward(values, fast)
    direct = cmplx(0.0_dp, 0.0_dp, dp)
    sign_transform = cmplx(0.0_dp, 0.0_dp, dp)
    permutation = [(permutation_index, permutation_index=1, degree)]
    do permutation_index = 1, order
      sign_transform = sign_transform+ &
           real(permutation_parity(permutation), dp)* &
           values(permutation_index)
      do block = 1, plan%number_of_blocks()
        dimension = plan%block_dimension(block)
        offset = plan%block_offset(block)
        call plan%representation_matrix(block, permutation, representation)
        direct(offset:offset+dimension**2-1) = &
             direct(offset:offset+dimension**2-1)+ &
             values(permutation_index)*reshape(representation, &
             [dimension**2])
      end do
      if (permutation_index < order) call next_permutation(permutation)
    end do

    call require_close(fast(1), sum(values), &
         'trivial-representation transform is incorrect')
    call require_close(fast(order), sign_transform, &
         'sign-representation transform is incorrect')

    scale = max(1.0_dp, maxval(abs(direct)))
    if (maxval(abs(fast-direct)) > 2.0e-11_dp*scale) then
      write(*, '(a,i0,2(1x,es24.16))') &
           'direct FFT mismatch at degree ', degree, &
           maxval(abs(fast-direct)), scale
      call fail('recursive symmetric-group transform is incorrect')
    end if

    input_norm = sum(abs(values)**2)
    fourier_norm = 0.0_dp
    do block = 1, plan%number_of_blocks()
      dimension = plan%block_dimension(block)
      offset = plan%block_offset(block)
      fourier_norm = fourier_norm+real(dimension, dp)* &
           sum(abs(fast(offset:offset+dimension**2-1))**2)
    end do
    fourier_norm = fourier_norm/real(order, dp)
    if (abs(fourier_norm-input_norm) > &
         2.0e-11_dp*max(1.0_dp, input_norm)) then
      write(*, '(a,i0,2(1x,es24.16))') &
           'Parseval mismatch at degree ', degree, fourier_norm, input_norm
      call fail('symmetric-group Parseval identity is violated')
    end if
  end subroutine check_direct_transform

  subroutine check_reusable_workspace(degree)
    integer, intent(in) :: degree

    type(symmetric_group_fft_t) :: plan
    type(symmetric_group_fft_workspace_t) :: workspace
    complex(dp), allocatable :: fresh(:), reusable(:), values(:)
    real(dp) :: scale
    integer :: order, value

    call plan%initialize(degree)
    call plan%initialize_workspace(workspace)
    order = plan%group_order()
    allocate(values(order), reusable(order))

    do value = 1, order
      values(value) = cmplx(sin(0.31_dp*real(value, dp)), &
           cos(0.29_dp*real(value, dp)), dp)
    end do
    call plan%forward_into(values, reusable, workspace)
    call plan%forward(values, fresh)
    scale = max(1.0_dp, maxval(abs(fresh)))
    if (maxval(abs(reusable-fresh)) > 2.0e-12_dp*scale) &
         call fail('preallocated symmetric-group transform is incorrect')

    ! A second, unrelated input checks that no state from the first transform
    ! survives when the same full-vector and block workspaces are reused.
    do value = 1, order
      values(value) = cmplx(cos(0.41_dp*real(value, dp))+0.01_dp*value, &
           sin(0.19_dp*real(value, dp))-0.02_dp*value, dp)
    end do
    call plan%forward_into(values, reusable, workspace)
    call plan%forward(values, fresh)
    scale = max(1.0_dp, maxval(abs(fresh)))
    if (maxval(abs(reusable-fresh)) > 2.0e-12_dp*scale) &
         call fail('reused symmetric-group FFT workspace retained state')
  end subroutine check_reusable_workspace

  subroutine check_parseval_transform(degree)
    integer, intent(in) :: degree

    type(symmetric_group_fft_t) :: plan
    complex(dp), allocatable :: coefficients(:), values(:)
    real(dp) :: fourier_norm, input_norm
    integer :: block, dimension, offset, order, value

    call plan%initialize(degree)
    order = plan%group_order()
    allocate(values(order))
    do value = 1, order
      values(value) = cmplx(sin(0.017_dp*real(value, dp)), &
           cos(0.011_dp*real(value, dp)), dp)
    end do
    call plan%forward(values, coefficients)

    input_norm = sum(abs(values)**2)
    fourier_norm = 0.0_dp
    do block = 1, plan%number_of_blocks()
      dimension = plan%block_dimension(block)
      offset = plan%block_offset(block)
      fourier_norm = fourier_norm+real(dimension, dp)* &
           sum(abs(coefficients(offset:offset+dimension**2-1))**2)
    end do
    fourier_norm = fourier_norm/real(order, dp)
    if (abs(fourier_norm-input_norm) > &
         5.0e-11_dp*max(1.0_dp, input_norm)) &
         call fail('large-degree symmetric-group Parseval mismatch')
  end subroutine check_parseval_transform

  subroutine check_selected_impulse(degree, impulse)
    integer, intent(in) :: degree, impulse

    type(symmetric_group_fft_t) :: plan
    complex(dp), allocatable :: coefficients(:), values(:)
    real(dp), allocatable :: representation(:, :)
    real(dp) :: scale
    integer :: block, dimension, offset, order, rank
    integer, allocatable :: permutation(:)

    call plan%initialize(degree)
    order = plan%group_order()
    call require(impulse >= 1 .and. impulse <= order, &
         'selected FFT impulse is out of range')
    allocate(values(order), permutation(degree))
    values = cmplx(0.0_dp, 0.0_dp, dp)
    values(impulse) = cmplx(1.0_dp, 0.0_dp, dp)
    call plan%forward(values, coefficients)

    permutation = [(rank, rank=1, degree)]
    do rank = 2, impulse
      call next_permutation(permutation)
    end do
    do block = 1, plan%number_of_blocks()
      dimension = plan%block_dimension(block)
      offset = plan%block_offset(block)
      call plan%representation_matrix(block, permutation, representation)
      scale = max(1.0_dp, maxval(abs(representation)))
      if (maxval(abs(coefficients(offset:offset+dimension**2-1)- &
           cmplx(reshape(representation, [dimension**2]), 0.0_dp, dp))) > &
           2.0e-11_dp*scale) &
           call fail('large-degree selected FFT impulse is incorrect')
    end do
  end subroutine check_selected_impulse

  subroutine check_representation_relations(degree)
    integer, intent(in) :: degree

    type(symmetric_group_fft_t) :: plan
    real(dp), allocatable :: first(:, :), second(:, :), third(:, :)
    real(dp), allocatable :: identity(:, :), left_product(:, :)
    real(dp), allocatable :: right_product(:, :)
    integer :: block, dimension, generator, row
    integer, allocatable :: permutation(:)

    call plan%initialize(degree)
    allocate(permutation(degree))
    do block = 1, plan%number_of_blocks()
      dimension = plan%block_dimension(block)
      allocate(identity(dimension, dimension))
      identity = 0.0_dp
      do row = 1, dimension
        identity(row, row) = 1.0_dp
      end do
      do generator = 1, degree-1
        call adjacent_permutation(degree, generator, permutation)
        call plan%representation_matrix(block, permutation, first)
        call require_matrix_close(matmul(first, first), identity, &
             'Young generator does not square to identity')
        call require_matrix_close(matmul(transpose(first), first), &
             identity, 'Young generator is not orthogonal')
      end do
      do generator = 1, degree-2
        call adjacent_permutation(degree, generator, permutation)
        call plan%representation_matrix(block, permutation, first)
        call adjacent_permutation(degree, generator+1, permutation)
        call plan%representation_matrix(block, permutation, second)
        left_product = matmul(first, matmul(second, first))
        right_product = matmul(second, matmul(first, second))
        call require_matrix_close(left_product, right_product, &
             'Young generators violate the braid relation')
      end do
      do generator = 1, degree-3
        call adjacent_permutation(degree, generator, permutation)
        call plan%representation_matrix(block, permutation, first)
        call adjacent_permutation(degree, generator+2, permutation)
        call plan%representation_matrix(block, permutation, third)
        call require_matrix_close(matmul(first, third), &
             matmul(third, first), &
             'distant Young generators do not commute')
      end do
      deallocate(identity)
    end do
  end subroutine check_representation_relations

  subroutine adjacent_permutation(degree, generator, permutation)
    integer, intent(in) :: degree, generator
    integer, intent(out) :: permutation(degree)

    integer :: temporary

    permutation = [(temporary, temporary=1, degree)]
    temporary = permutation(generator)
    permutation(generator) = permutation(generator+1)
    permutation(generator+1) = temporary
  end subroutine adjacent_permutation

  integer function permutation_parity(permutation) result(parity)
    integer, intent(in) :: permutation(:)

    integer :: inversions, left, right

    inversions = 0
    do left = 1, size(permutation)-1
      do right = left+1, size(permutation)
        if (permutation(left) > permutation(right)) &
             inversions = inversions+1
      end do
    end do
    parity = merge(1, -1, mod(inversions, 2) == 0)
  end function permutation_parity

  subroutine next_permutation(permutation)
    integer, intent(inout) :: permutation(:)

    integer :: left, right, temporary

    left = size(permutation)-1
    do while (left >= 1)
      if (permutation(left) < permutation(left+1)) exit
      left = left-1
    end do
    if (left < 1) call fail('advanced beyond the final permutation')
    right = size(permutation)
    do while (permutation(right) < permutation(left))
      right = right-1
    end do
    temporary = permutation(left)
    permutation(left) = permutation(right)
    permutation(right) = temporary
    permutation(left+1:) = permutation(size(permutation):left+1:-1)
  end subroutine next_permutation

  integer function factorial(number) result(value)
    integer, intent(in) :: number

    integer :: factor

    value = 1
    do factor = 2, number
      value = value*factor
    end do
  end function factorial

  subroutine require_matrix_close(actual, expected, message)
    real(dp), intent(in) :: actual(:, :), expected(:, :)
    character(len=*), intent(in) :: message

    if (maxval(abs(actual-expected)) > 5.0e-12_dp) call fail(message)
  end subroutine require_matrix_close

  subroutine require_close(actual, expected, message)
    complex(dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: message

    if (abs(actual-expected) > 5.0e-12_dp) call fail(message)
  end subroutine require_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) call fail(message)
  end subroutine require

end program test_symmetric_group_fft
