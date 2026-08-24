module symmetric_group_fft
  use ampligluon_common, only: dp, fail, i64
  implicit none
  private

  ! The forward transform is unnormalised:
  !
  !   F_lambda = sum_sigma f(sigma) rho_lambda(sigma).
  !
  ! Function values use lexicographic image-permutation order.  Fourier
  ! matrices are flattened column-major into contiguous partition blocks.

  type :: young_irrep_t
    integer :: dimension = 0
    integer :: coefficient_offset = 0
    integer, allocatable :: shape(:)
    integer, allocatable :: child_index(:)
    integer, allocatable :: child_basis_offset(:)
    integer, allocatable :: tableau_row(:, :)
    integer, allocatable :: tableau_column(:, :)
    integer, allocatable :: generator_partner(:, :)
    real(dp), allocatable :: generator_diagonal(:, :)
    real(dp), allocatable :: generator_mixing(:, :)
  end type young_irrep_t

  type :: symmetric_group_level_t
    integer :: degree = 0
    integer :: group_order = 0
    integer :: maximum_dimension = 0
    integer, allocatable :: recursive_to_lexicographic(:)
    type(young_irrep_t), allocatable :: irrep(:)
  end type symmetric_group_level_t

  ! Reusable storage for forward_into.  Keeping this separate from the plan
  ! lets one immutable plan serve callers that need independent workspaces.
  type, public :: symmetric_group_fft_workspace_t
    private
    logical :: initialized = .false.
    integer :: degree_value = 0
    integer :: group_order_value = 0
    complex(dp), allocatable :: buffer(:)
    complex(dp), allocatable :: block(:, :)
  end type symmetric_group_fft_workspace_t

  type, public :: symmetric_group_fft_t
    private
    logical :: initialized = .false.
    integer :: degree_value = 0
    type(symmetric_group_level_t), allocatable :: level(:)
  contains
    procedure, public :: block_dimension
    procedure, public :: block_layout
    procedure, public :: block_offset
    procedure, public :: block_shape
    procedure, public :: degree
    procedure, public :: forward
    procedure, public :: forward_into
    procedure, public :: group_order
    procedure, public :: initialize
    procedure, public :: initialize_workspace
    procedure, public :: number_of_blocks
    procedure, public :: representation_matrix
  end type symmetric_group_fft_t

contains

  subroutine initialize(this, degree)
    class(symmetric_group_fft_t), intent(inout) :: this
    integer, intent(in) :: degree

    integer :: current_degree
    integer(i64) :: order

    this%initialized = .false.
    this%degree_value = 0
    if (allocated(this%level)) deallocate(this%level)
    if (degree < 1) call fail('the symmetric-group degree must be positive')

    allocate(this%level(degree))
    order = 1_i64
    do current_degree = 1, degree
      if (order > int(huge(0), i64)) &
           call fail('the symmetric-group order exceeds array-index capacity')
      call initialize_level(this%level, current_degree, int(order))
      if (current_degree > 1) &
           call connect_branching(this%level(current_degree), &
           this%level(current_degree-1))
      call initialize_recursive_order(this%level, current_degree)
      if (current_degree < degree) then
        if (order > huge(order)/int(current_degree+1, i64)) &
             call fail('the symmetric-group order overflows a 64-bit integer')
        order = order*int(current_degree+1, i64)
      end if
    end do

    this%degree_value = degree
    this%initialized = .true.
  end subroutine initialize

  subroutine initialize_level(levels, degree, order)
    type(symmetric_group_level_t), intent(inout) :: levels(:)
    integer, intent(in) :: degree, order

    integer :: block, coefficient_offset, number_of_partitions
    integer, allocatable :: partitions(:, :)

    number_of_partitions = count_partitions(degree)
    allocate(partitions(degree, number_of_partitions))
    call enumerate_partitions(degree, partitions)

    levels(degree)%degree = degree
    levels(degree)%group_order = order
    allocate(levels(degree)%irrep(number_of_partitions))
    do block = 1, number_of_partitions
      allocate(levels(degree)%irrep(block)%shape(degree))
      levels(degree)%irrep(block)%shape = partitions(:, block)
    end do

    if (degree == 1) then
      levels(1)%irrep(1)%dimension = 1
      levels(1)%irrep(1)%coefficient_offset = 1
      allocate(levels(1)%irrep(1)%child_index(0))
      allocate(levels(1)%irrep(1)%child_basis_offset(0))
      allocate(levels(1)%irrep(1)%tableau_row(1, 1))
      allocate(levels(1)%irrep(1)%tableau_column(1, 1))
      allocate(levels(1)%irrep(1)%generator_partner(0, 1))
      allocate(levels(1)%irrep(1)%generator_diagonal(0, 1))
      allocate(levels(1)%irrep(1)%generator_mixing(0, 1))
      levels(1)%irrep(1)%tableau_row = 1
      levels(1)%irrep(1)%tableau_column = 1
      levels(1)%maximum_dimension = 1
      return
    end if

    coefficient_offset = 1
    do block = 1, number_of_partitions
      call initialize_irrep(levels(degree)%irrep(block), &
           levels(degree-1), degree)
      levels(degree)%irrep(block)%coefficient_offset = coefficient_offset
      coefficient_offset = coefficient_offset+ &
           levels(degree)%irrep(block)%dimension**2
      levels(degree)%maximum_dimension = max( &
           levels(degree)%maximum_dimension, &
           levels(degree)%irrep(block)%dimension)
    end do
    if (coefficient_offset-1 /= order) &
         call fail('Young dimensions do not square to the group order')
  end subroutine initialize_level

  subroutine initialize_irrep(irrep, child_level, degree)
    type(young_irrep_t), intent(inout) :: irrep
    type(symmetric_group_level_t), intent(in) :: child_level
    integer, intent(in) :: degree

    integer :: basis, branch, child, child_dimension, corner, dimension
    integer :: number_of_corners

    number_of_corners = 0
    do corner = 1, degree
      if (is_removable_corner(irrep%shape, corner)) &
           number_of_corners = number_of_corners+1
    end do
    allocate(irrep%child_index(number_of_corners))
    allocate(irrep%child_basis_offset(number_of_corners))

    dimension = 0
    branch = 0
    do corner = 1, degree
      if (.not. is_removable_corner(irrep%shape, corner)) cycle
      branch = branch+1
      child = find_child_irrep(irrep%shape, corner, child_level)
      irrep%child_index(branch) = child
      irrep%child_basis_offset(branch) = dimension+1
      dimension = dimension+child_level%irrep(child)%dimension
    end do
    irrep%dimension = dimension

    allocate(irrep%tableau_row(degree, dimension))
    allocate(irrep%tableau_column(degree, dimension))
    branch = 0
    do corner = 1, degree
      if (.not. is_removable_corner(irrep%shape, corner)) cycle
      branch = branch+1
      child = irrep%child_index(branch)
      child_dimension = child_level%irrep(child)%dimension
      basis = irrep%child_basis_offset(branch)
      irrep%tableau_row(1:degree-1, basis:basis+child_dimension-1) = &
           child_level%irrep(child)%tableau_row
      irrep%tableau_column(1:degree-1, &
           basis:basis+child_dimension-1) = &
           child_level%irrep(child)%tableau_column
      irrep%tableau_row(degree, basis:basis+child_dimension-1) = corner
      irrep%tableau_column(degree, basis:basis+child_dimension-1) = &
           irrep%shape(corner)
    end do

    call initialize_generators(irrep, degree)
  end subroutine initialize_irrep

  subroutine initialize_generators(irrep, degree)
    type(young_irrep_t), intent(inout) :: irrep
    integer, intent(in) :: degree

    integer :: axial_distance, basis, generator, partner
    integer(i64) :: target_code
    integer(i64), allocatable :: sorted_code(:), tableau_code(:)
    integer, allocatable :: sorted_index(:)

    allocate(irrep%generator_partner(degree-1, irrep%dimension))
    allocate(irrep%generator_diagonal(degree-1, irrep%dimension))
    allocate(irrep%generator_mixing(degree-1, irrep%dimension))
    allocate(tableau_code(irrep%dimension))
    allocate(sorted_code(irrep%dimension), sorted_index(irrep%dimension))
    do basis = 1, irrep%dimension
      tableau_code(basis) = encode_tableau(irrep%tableau_row(:, basis))
      sorted_code(basis) = tableau_code(basis)
      sorted_index(basis) = basis
    end do
    call sort_codes(sorted_code, sorted_index, 1, irrep%dimension)

    do basis = 1, irrep%dimension
      do generator = 1, degree-1
        axial_distance = &
             irrep%tableau_column(generator+1, basis)- &
             irrep%tableau_row(generator+1, basis)- &
             irrep%tableau_column(generator, basis)+ &
             irrep%tableau_row(generator, basis)
        if (axial_distance == 0) &
             call fail('zero axial distance in a standard Young tableau')
        irrep%generator_diagonal(generator, basis) = &
             1.0_dp/real(axial_distance, dp)
        if (abs(axial_distance) == 1) then
          partner = basis
        else
          target_code = swapped_tableau_code( &
               irrep%tableau_row(:, basis), generator)
          partner = lookup_code(sorted_code, sorted_index, target_code)
          if (partner == 0) &
               call fail('adjacent tableau partner was not generated')
        end if
        irrep%generator_partner(generator, basis) = partner
        irrep%generator_mixing(generator, basis) = sqrt(max(0.0_dp, &
             1.0_dp-irrep%generator_diagonal(generator, basis)**2))
      end do
    end do
  end subroutine initialize_generators

  subroutine connect_branching(level, child_level)
    type(symmetric_group_level_t), intent(in) :: level, child_level

    integer :: block, branch, child, child_dimension

    do block = 1, size(level%irrep)
      do branch = 1, size(level%irrep(block)%child_index)
        child = level%irrep(block)%child_index(branch)
        child_dimension = child_level%irrep(child)%dimension
        if (level%irrep(block)%child_basis_offset(branch)+ &
             child_dimension-1 > level%irrep(block)%dimension) &
             call fail('invalid Young branching offset')
      end do
    end do
  end subroutine connect_branching

  subroutine initialize_recursive_order(levels, degree)
    type(symmetric_group_level_t), intent(inout) :: levels(:)
    integer, intent(in) :: degree

    integer :: child_order, coset, recursive_index
    integer, allocatable :: child_permutation(:), permutation(:)

    allocate(levels(degree)%recursive_to_lexicographic( &
         levels(degree)%group_order))
    if (degree == 1) then
      levels(1)%recursive_to_lexicographic = 1
      return
    end if

    child_order = levels(degree-1)%group_order
    allocate(child_permutation(degree-1), permutation(degree))
    do coset = 1, degree
      do recursive_index = 1, child_order
        call lexicographic_permutation(degree-1, &
             levels(degree-1)%recursive_to_lexicographic( &
             recursive_index), child_permutation)
        where (child_permutation < coset)
          permutation(1:degree-1) = child_permutation
        elsewhere
          permutation(1:degree-1) = child_permutation+1
        end where
        permutation(degree) = coset
        levels(degree)%recursive_to_lexicographic( &
             (coset-1)*child_order+recursive_index) = &
             lexicographic_rank(permutation)
      end do
    end do
  end subroutine initialize_recursive_order

  subroutine forward(this, function_values, coefficients)
    class(symmetric_group_fft_t), intent(in) :: this
    complex(dp), intent(in) :: function_values(:)
    complex(dp), allocatable, intent(out) :: coefficients(:)

    type(symmetric_group_fft_workspace_t) :: workspace
    integer :: order

    call require_initialized(this)
    order = this%level(this%degree_value)%group_order
    if (size(function_values) /= order) &
         call fail('wrong number of symmetric-group function values')

    allocate(coefficients(order))
    call this%initialize_workspace(workspace)
    call this%forward_into(function_values, coefficients, workspace)
  end subroutine forward

  subroutine initialize_workspace(this, workspace)
    class(symmetric_group_fft_t), intent(in) :: this
    type(symmetric_group_fft_workspace_t), intent(inout) :: workspace

    integer :: current_degree, maximum_dimension, order

    call require_initialized(this)
    workspace%initialized = .false.
    workspace%degree_value = 0
    workspace%group_order_value = 0
    if (allocated(workspace%buffer)) deallocate(workspace%buffer)
    if (allocated(workspace%block)) deallocate(workspace%block)

    order = this%level(this%degree_value)%group_order
    maximum_dimension = 1
    do current_degree = 1, this%degree_value
      maximum_dimension = max(maximum_dimension, &
           this%level(current_degree)%maximum_dimension)
    end do
    allocate(workspace%buffer(order))
    allocate(workspace%block(maximum_dimension, maximum_dimension))
    workspace%degree_value = this%degree_value
    workspace%group_order_value = order
    workspace%initialized = .true.
  end subroutine initialize_workspace

  subroutine forward_into(this, function_values, coefficients, workspace)
    class(symmetric_group_fft_t), intent(in) :: this
    complex(dp), intent(in) :: function_values(:)
    complex(dp), intent(out) :: coefficients(:)
    type(symmetric_group_fft_workspace_t), intent(inout) :: workspace

    integer :: current_degree, group, group_start, order, value
    logical :: input_is_workspace

    call require_initialized(this)
    order = this%level(this%degree_value)%group_order
    if (size(function_values) /= order) &
         call fail('wrong number of symmetric-group function values')
    if (size(coefficients) /= order) &
         call fail('wrong number of symmetric-group Fourier coefficients')
    call require_compatible_workspace(this, workspace)

    ! A one-time reorder exposes the nested left-coset decomposition
    ! S_m = union_i (s_i ... s_{m-1}) S_{m-1} at every subgroup level.
    do value = 1, order
      workspace%buffer(value) = function_values( &
           this%level(this%degree_value)%recursive_to_lexicographic(value))
    end do
    input_is_workspace = .true.
    do current_degree = 2, this%degree_value
      do group = 1, order/this%level(current_degree)%group_order
        group_start = (group-1)*this%level(current_degree)%group_order+1
        if (input_is_workspace) then
          call transform_group(this%level(current_degree), &
               this%level(current_degree-1), &
               workspace%buffer(group_start:group_start+ &
               this%level(current_degree)%group_order-1), &
               coefficients(group_start:group_start+ &
               this%level(current_degree)%group_order-1), workspace%block)
        else
          call transform_group(this%level(current_degree), &
               this%level(current_degree-1), &
               coefficients(group_start:group_start+ &
               this%level(current_degree)%group_order-1), &
               workspace%buffer(group_start:group_start+ &
               this%level(current_degree)%group_order-1), workspace%block)
        end if
      end do
      input_is_workspace = .not. input_is_workspace
    end do

    if (input_is_workspace) then
      call transpose_fourier_blocks(this%level(this%degree_value), &
           workspace%buffer, coefficients)
    else
      call transpose_fourier_blocks_in_place( &
           this%level(this%degree_value), coefficients)
    end if
  end subroutine forward_into

  subroutine transform_group(level, child_level, child_transforms, &
       transform, workspace)
    type(symmetric_group_level_t), intent(in) :: level, child_level
    complex(dp), intent(in) :: child_transforms(:)
    complex(dp), intent(out) :: transform(:)
    complex(dp), intent(inout) :: workspace(:, :)

    integer :: basis_offset, block, branch, child, child_dimension
    integer :: child_offset, column, coset, dimension, generator, row
    integer :: transform_offset

    transform = cmplx(0.0_dp, 0.0_dp, dp)
    do coset = 1, level%degree
      child_offset = (coset-1)*child_level%group_order
      do block = 1, size(level%irrep)
        dimension = level%irrep(block)%dimension
        workspace(1:dimension, 1:dimension) = &
             cmplx(0.0_dp, 0.0_dp, dp)
        do branch = 1, size(level%irrep(block)%child_index)
          ! Restriction to S_(m-1) is block diagonal because tableaux in
          ! each parent basis are grouped by the shape left after removing m.
          child = level%irrep(block)%child_index(branch)
          child_dimension = child_level%irrep(child)%dimension
          basis_offset = level%irrep(block)%child_basis_offset(branch)
          do column = 1, child_dimension
            do row = 1, child_dimension
              workspace(basis_offset+row-1, basis_offset+column-1) = &
                   child_transforms(child_offset+ &
                   child_level%irrep(child)%coefficient_offset+ &
                   (column-1)*child_dimension+row-1)
            end do
          end do
        end do
        ! Left multiplication in this order produces
        ! rho(s_i)...rho(s_(m-1)) times the restricted child transform.
        do generator = level%degree-1, coset, -1
          call right_apply_generator_transposed( &
               level%irrep(block), generator, &
               workspace(1:dimension, 1:dimension))
        end do
        transform_offset = level%irrep(block)%coefficient_offset
        do column = 1, dimension
          do row = 1, dimension
            transform(transform_offset+(column-1)*dimension+row-1) = &
                 transform(transform_offset+(column-1)*dimension+row-1)+ &
                 workspace(row, column)
          end do
        end do
      end do
    end do
  end subroutine transform_group

  subroutine right_apply_generator_transposed(irrep, generator, matrix)
    type(young_irrep_t), intent(in) :: irrep
    integer, intent(in) :: generator
    complex(dp), intent(inout) :: matrix(:, :)

    complex(dp) :: first_value, second_value
    real(dp) :: first_diagonal, first_imaginary, first_real, mixing
    real(dp) :: second_diagonal, second_imaginary, second_real
    integer :: first, row, second

    ! Intermediate Fourier blocks are stored transposed.  Left multiplication
    ! by a real symmetric Young generator therefore becomes right
    ! multiplication here, making the hot row traversal unit-stride.
    do first = 1, irrep%dimension
      second = irrep%generator_partner(generator, first)
      first_diagonal = irrep%generator_diagonal(generator, first)
      if (second == first) then
        if (first_diagonal > 0.0_dp) cycle
        do row = 1, irrep%dimension
          first_value = matrix(row, first)
          matrix(row, first) = cmplx(-real(first_value, dp), &
               -aimag(first_value), dp)
        end do
      else if (second > first) then
        mixing = irrep%generator_mixing(generator, first)
        second_diagonal = -first_diagonal
        do row = 1, irrep%dimension
          first_value = matrix(row, first)
          second_value = matrix(row, second)
          first_real = real(first_value, dp)
          first_imaginary = aimag(first_value)
          second_real = real(second_value, dp)
          second_imaginary = aimag(second_value)
          matrix(row, first) = cmplx(first_diagonal*first_real+ &
               mixing*second_real, first_diagonal*first_imaginary+ &
               mixing*second_imaginary, dp)
          matrix(row, second) = cmplx(mixing*first_real+ &
               second_diagonal*second_real, mixing*first_imaginary+ &
               second_diagonal*second_imaginary, dp)
        end do
      end if
    end do
  end subroutine right_apply_generator_transposed

  subroutine transpose_fourier_blocks(level, transposed, coefficients)
    type(symmetric_group_level_t), intent(in) :: level
    complex(dp), intent(in) :: transposed(:)
    complex(dp), intent(out) :: coefficients(:)

    integer :: block, column, dimension, offset, row

    do block = 1, size(level%irrep)
      dimension = level%irrep(block)%dimension
      offset = level%irrep(block)%coefficient_offset
      do column = 1, dimension
        do row = 1, dimension
          coefficients(offset+(column-1)*dimension+row-1) = &
               transposed(offset+(row-1)*dimension+column-1)
        end do
      end do
    end do
  end subroutine transpose_fourier_blocks

  subroutine transpose_fourier_blocks_in_place(level, coefficients)
    type(symmetric_group_level_t), intent(in) :: level
    complex(dp), intent(inout) :: coefficients(:)

    complex(dp) :: temporary
    integer :: block, column, dimension, first, offset, row, second

    do block = 1, size(level%irrep)
      dimension = level%irrep(block)%dimension
      offset = level%irrep(block)%coefficient_offset
      do column = 1, dimension
        do row = column+1, dimension
          first = offset+(column-1)*dimension+row-1
          second = offset+(row-1)*dimension+column-1
          temporary = coefficients(first)
          coefficients(first) = coefficients(second)
          coefficients(second) = temporary
        end do
      end do
    end do
  end subroutine transpose_fourier_blocks_in_place

  subroutine left_apply_generator(irrep, generator, matrix)
    type(young_irrep_t), intent(in) :: irrep
    integer, intent(in) :: generator
    complex(dp), intent(inout) :: matrix(:, :)

    complex(dp) :: first_value, second_value
    real(dp) :: mixing
    integer :: column, first, second

    do first = 1, irrep%dimension
      second = irrep%generator_partner(generator, first)
      if (second == first) then
        do column = 1, irrep%dimension
          matrix(first, column) = &
               irrep%generator_diagonal(generator, first)* &
               matrix(first, column)
        end do
      else if (second > first) then
        mixing = irrep%generator_mixing(generator, first)
        do column = 1, irrep%dimension
          first_value = matrix(first, column)
          second_value = matrix(second, column)
          matrix(first, column) = &
               irrep%generator_diagonal(generator, first)*first_value+ &
               mixing*second_value
          matrix(second, column) = mixing*first_value+ &
               irrep%generator_diagonal(generator, second)*second_value
        end do
      end if
    end do
  end subroutine left_apply_generator

  subroutine representation_matrix(this, block, permutation, matrix)
    class(symmetric_group_fft_t), intent(in) :: this
    integer, intent(in) :: block, permutation(:)
    real(dp), allocatable, intent(out) :: matrix(:, :)

    complex(dp), allocatable :: complex_matrix(:, :)
    integer :: dimension, left, pass, temporary
    integer, allocatable :: work_permutation(:)
    logical :: changed

    call check_block(this, block)
    call validate_permutation(permutation, this%degree_value)
    dimension = this%level(this%degree_value)%irrep(block)%dimension
    allocate(complex_matrix(dimension, dimension))
    complex_matrix = cmplx(0.0_dp, 0.0_dp, dp)
    do left = 1, dimension
      complex_matrix(left, left) = cmplx(1.0_dp, 0.0_dp, dp)
    end do

    work_permutation = permutation
    do pass = 1, this%degree_value-1
      changed = .false.
      do left = 1, this%degree_value-pass
        if (work_permutation(left) <= work_permutation(left+1)) cycle
        temporary = work_permutation(left)
        work_permutation(left) = work_permutation(left+1)
        work_permutation(left+1) = temporary
        call left_apply_generator( &
             this%level(this%degree_value)%irrep(block), left, &
             complex_matrix)
        changed = .true.
      end do
      if (.not. changed) exit
    end do
    if (any(work_permutation /= [(left, left=1, this%degree_value)])) &
         call fail('failed to reduce a permutation to adjacent generators')

    allocate(matrix(dimension, dimension))
    matrix = real(complex_matrix, dp)
  end subroutine representation_matrix

  integer function degree(this) result(value)
    class(symmetric_group_fft_t), intent(in) :: this

    call require_initialized(this)
    value = this%degree_value
  end function degree

  integer function group_order(this) result(value)
    class(symmetric_group_fft_t), intent(in) :: this

    call require_initialized(this)
    value = this%level(this%degree_value)%group_order
  end function group_order

  integer function number_of_blocks(this) result(number)
    class(symmetric_group_fft_t), intent(in) :: this

    call require_initialized(this)
    number = size(this%level(this%degree_value)%irrep)
  end function number_of_blocks

  integer function block_dimension(this, block) result(dimension)
    class(symmetric_group_fft_t), intent(in) :: this
    integer, intent(in) :: block

    call check_block(this, block)
    dimension = this%level(this%degree_value)%irrep(block)%dimension
  end function block_dimension

  integer function block_offset(this, block) result(offset)
    class(symmetric_group_fft_t), intent(in) :: this
    integer, intent(in) :: block

    call check_block(this, block)
    offset = this%level(this%degree_value)%irrep(block)%coefficient_offset
  end function block_offset

  function block_shape(this, block) result(shape)
    class(symmetric_group_fft_t), intent(in) :: this
    integer, intent(in) :: block
    integer, allocatable :: shape(:)

    call check_block(this, block)
    shape = this%level(this%degree_value)%irrep(block)%shape
  end function block_shape

  subroutine block_layout(this, offsets, dimensions)
    class(symmetric_group_fft_t), intent(in) :: this
    integer, allocatable, intent(out) :: offsets(:), dimensions(:)

    integer :: block

    call require_initialized(this)
    allocate(offsets(this%number_of_blocks()))
    allocate(dimensions(this%number_of_blocks()))
    do block = 1, this%number_of_blocks()
      offsets(block) = this%block_offset(block)
      dimensions(block) = this%block_dimension(block)
    end do
  end subroutine block_layout

  subroutine require_initialized(this)
    class(symmetric_group_fft_t), intent(in) :: this

    if (.not. this%initialized) &
         call fail('the symmetric-group FFT is not initialized')
  end subroutine require_initialized

  subroutine require_compatible_workspace(this, workspace)
    class(symmetric_group_fft_t), intent(in) :: this
    type(symmetric_group_fft_workspace_t), intent(in) :: workspace

    integer :: current_degree, maximum_dimension, order

    if (.not. workspace%initialized) &
         call fail('the symmetric-group FFT workspace is not initialized')
    order = this%level(this%degree_value)%group_order
    if (workspace%degree_value /= this%degree_value .or. &
         workspace%group_order_value /= order) &
         call fail('the symmetric-group FFT workspace belongs to another plan')
    if (.not. allocated(workspace%buffer) .or. &
         size(workspace%buffer) /= order) &
         call fail('the symmetric-group FFT workspace buffer has the wrong size')
    if (.not. allocated(workspace%block)) &
         call fail('the symmetric-group FFT block workspace is not allocated')
    maximum_dimension = 1
    do current_degree = 1, this%degree_value
      maximum_dimension = max(maximum_dimension, &
           this%level(current_degree)%maximum_dimension)
    end do
    if (size(workspace%block, 1) < maximum_dimension .or. &
         size(workspace%block, 2) < maximum_dimension) &
         call fail('the symmetric-group FFT block workspace is too small')
  end subroutine require_compatible_workspace

  subroutine check_block(this, block)
    class(symmetric_group_fft_t), intent(in) :: this
    integer, intent(in) :: block

    call require_initialized(this)
    if (block < 1 .or. block > &
         size(this%level(this%degree_value)%irrep)) &
         call fail('symmetric-group Fourier block index is out of range')
  end subroutine check_block

  integer function count_partitions(number) result(count)
    integer, intent(in) :: number

    integer :: largest

    count = 0
    do largest = 1, number
      count = count+count_partitions_with_maximum(number-largest, largest)
    end do
  end function count_partitions

  recursive integer function count_partitions_with_maximum( &
       remaining, maximum) result(count)
    integer, intent(in) :: remaining, maximum

    integer :: part

    if (remaining == 0) then
      count = 1
      return
    end if
    count = 0
    do part = min(remaining, maximum), 1, -1
      count = count+count_partitions_with_maximum(remaining-part, part)
    end do
  end function count_partitions_with_maximum

  subroutine enumerate_partitions(number, partitions)
    integer, intent(in) :: number
    integer, intent(out) :: partitions(:, :)

    integer :: cursor
    integer, allocatable :: current(:)

    allocate(current(number))
    current = 0
    partitions = 0
    cursor = 0
    call append_partitions(number, number, 1, current, partitions, cursor)
    if (cursor /= size(partitions, 2)) &
         call fail('partition enumeration count mismatch')
  end subroutine enumerate_partitions

  recursive subroutine append_partitions(remaining, maximum, position, &
       current, partitions, cursor)
    integer, intent(in) :: remaining, maximum, position
    integer, intent(inout) :: current(:), partitions(:, :)
    integer, intent(inout) :: cursor

    integer :: part

    if (remaining == 0) then
      cursor = cursor+1
      partitions(:, cursor) = current
      return
    end if
    do part = min(remaining, maximum), 1, -1
      current(position) = part
      call append_partitions(remaining-part, part, position+1, current, &
           partitions, cursor)
      current(position) = 0
    end do
  end subroutine append_partitions

  logical function is_removable_corner(shape, row) result(removable)
    integer, intent(in) :: shape(:), row

    if (row < size(shape)) then
      removable = shape(row) > shape(row+1)
    else
      removable = shape(row) > 0
    end if
  end function is_removable_corner

  integer function find_child_irrep(shape, corner, child_level) &
       result(child)
    integer, intent(in) :: shape(:), corner
    type(symmetric_group_level_t), intent(in) :: child_level

    integer :: candidate(size(shape)-1), block

    candidate = shape(1:size(shape)-1)
    if (corner < size(shape)) candidate(corner) = candidate(corner)-1
    do block = 1, size(child_level%irrep)
      if (all(candidate == child_level%irrep(block)%shape)) then
        child = block
        return
      end if
    end do
    call fail('Young branching child shape was not found')
    child = 0
  end function find_child_irrep

  integer(i64) function encode_tableau(rows) result(code)
    integer, intent(in) :: rows(:)

    integer :: entry
    integer(i64) :: base

    base = int(size(rows)+1, i64)
    code = 0_i64
    do entry = 1, size(rows)
      if (code > (huge(code)-int(rows(entry), i64))/base) &
           call fail('Young-tableau code overflows a 64-bit integer')
      code = code*base+int(rows(entry), i64)
    end do
  end function encode_tableau

  integer(i64) function swapped_tableau_code(rows, generator) result(code)
    integer, intent(in) :: rows(:), generator

    integer :: entry, row
    integer(i64) :: base

    base = int(size(rows)+1, i64)
    code = 0_i64
    do entry = 1, size(rows)
      if (entry == generator) then
        row = rows(generator+1)
      else if (entry == generator+1) then
        row = rows(generator)
      else
        row = rows(entry)
      end if
      if (code > (huge(code)-int(row, i64))/base) &
           call fail('Young-tableau code overflows a 64-bit integer')
      code = code*base+int(row, i64)
    end do
  end function swapped_tableau_code

  recursive subroutine sort_codes(codes, indices, first, last)
    integer(i64), intent(inout) :: codes(:)
    integer, intent(inout) :: indices(:)
    integer, intent(in) :: first, last

    integer :: left, right, temporary_index
    integer(i64) :: pivot, temporary_code

    if (first >= last) return
    left = first
    right = last
    pivot = codes((first+last)/2)
    do
      do while (left <= last)
        if (codes(left) >= pivot) exit
        left = left+1
      end do
      do while (right >= first)
        if (codes(right) <= pivot) exit
        right = right-1
      end do
      if (left > right) exit
      temporary_code = codes(left)
      codes(left) = codes(right)
      codes(right) = temporary_code
      temporary_index = indices(left)
      indices(left) = indices(right)
      indices(right) = temporary_index
      left = left+1
      right = right-1
    end do
    if (first < right) call sort_codes(codes, indices, first, right)
    if (left < last) call sort_codes(codes, indices, left, last)
  end subroutine sort_codes

  integer function lookup_code(codes, indices, target) result(index)
    integer(i64), intent(in) :: codes(:), target
    integer, intent(in) :: indices(:)

    integer :: middle, lower, upper

    lower = 1
    upper = size(codes)
    do while (lower <= upper)
      middle = (lower+upper)/2
      if (codes(middle) < target) then
        lower = middle+1
      else if (codes(middle) > target) then
        upper = middle-1
      else
        index = indices(middle)
        return
      end if
    end do
    index = 0
  end function lookup_code

  subroutine lexicographic_permutation(degree, rank, permutation)
    integer, intent(in) :: degree, rank
    integer, intent(out) :: permutation(degree)

    integer :: available(degree), available_count, digit, factorial
    integer :: position, remainder

    if (rank < 1) call fail('lexicographic permutation rank is invalid')
    available = [(position, position=1, degree)]
    available_count = degree
    remainder = rank-1
    do position = 1, degree
      factorial = factorial_integer(degree-position)
      digit = remainder/factorial+1
      remainder = modulo(remainder, factorial)
      if (digit > available_count) &
           call fail('lexicographic permutation rank is too large')
      permutation(position) = available(digit)
      available(digit:available_count-1) = available(digit+1:available_count)
      available_count = available_count-1
    end do
  end subroutine lexicographic_permutation

  integer function lexicographic_rank(permutation) result(rank)
    integer, intent(in) :: permutation(:)

    integer :: left, right, smaller

    rank = 1
    do left = 1, size(permutation)-1
      smaller = 0
      do right = left+1, size(permutation)
        if (permutation(right) < permutation(left)) smaller = smaller+1
      end do
      rank = rank+smaller*factorial_integer(size(permutation)-left)
    end do
  end function lexicographic_rank

  integer function factorial_integer(number) result(value)
    integer, intent(in) :: number

    integer :: factor

    value = 1
    do factor = 2, number
      if (value > huge(value)/factor) &
           call fail('factorial exceeds array-index capacity')
      value = value*factor
    end do
  end function factorial_integer

  subroutine validate_permutation(permutation, degree)
    integer, intent(in) :: permutation(:), degree

    logical :: seen(degree)
    integer :: position

    if (size(permutation) /= degree) &
         call fail('wrong permutation degree for the Fourier plan')
    seen = .false.
    do position = 1, degree
      if (permutation(position) < 1 .or. permutation(position) > degree) &
           call fail('permutation entry is out of range')
      if (seen(permutation(position))) &
           call fail('duplicate entry in a permutation')
      seen(permutation(position)) = .true.
    end do
  end subroutine validate_permutation

end module symmetric_group_fft
