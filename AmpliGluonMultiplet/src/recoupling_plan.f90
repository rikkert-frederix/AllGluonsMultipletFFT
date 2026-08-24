module recoupling_plan
  use iso_c_binding, only: c_bool
  use ampligluon_multiplet_kinds, only: dp, fail, i8
  use multiplet_paths, only: conjugate_label, path_catalog_t, path_space_t
  use wigner_table, only: wigner_table_t
  implicit none
  private

  real(dp), parameter :: sparse_zero = 1.0e-15_dp

  type, public :: sparse_matrix_t
    integer :: size = 0
    integer :: number_of_values = 0
    logical :: is_signed_permutation = .false.
    integer, allocatable :: column_start(:)
    integer, allocatable :: row(:)
    real(dp), allocatable :: value(:)
    integer, allocatable :: signed_permutation(:)
  end type sparse_matrix_t

  type, public :: swap_set_t
    type(sparse_matrix_t), allocatable :: matrices(:)
  end type swap_set_t

  type, public :: path_support_t
    integer, allocatable :: active_path(:)
  end type path_support_t

  type, public :: ordering_node_t
    integer, allocatable :: order(:)
    integer, allocatable :: active_path(:)
    integer :: next_node = 0
    integer :: swap_position = 0
    integer :: inversions = 0
    integer :: layer_slot = 0
    logical :: destination_preinitialized = .false.
  end type ordering_node_t

  type, public :: partition_t
    integer(i8) :: left_pattern = 0_i8
    integer :: left_size = 0
    integer :: right_size = 0
    integer :: initial_node = 0
  end type partition_t

  type, public :: seed_map_t
    integer :: left_size = 0
    integer :: right_size = 0
    integer, allocatable :: signed_parent_path(:, :)
  end type seed_map_t

  type, public :: recoupling_plan_t
    integer :: length = 0
    integer :: number_of_paths = 0
    integer :: number_of_nodes = 0
    integer :: number_of_partitions = 0
    integer :: canonical_node = 0
    integer :: max_inversions = 0
    integer :: max_layer_width = 0
    type(ordering_node_t), allocatable :: nodes(:)
    type(partition_t), allocatable :: partitions(:)
    type(seed_map_t), allocatable :: seeds(:)
    integer, allocatable :: layer_node_start(:)
    integer, allocatable :: layer_nodes(:)
    integer, allocatable :: layer_partition_start(:)
    integer, allocatable :: layer_partitions(:)
  end type recoupling_plan_t

  type, public :: recoupling_system_t
    integer :: max_length = 0
    logical :: swaps_built = .false.
    logical :: plans_built = .false.
    type(swap_set_t), allocatable :: swaps(:)
    type(recoupling_plan_t), allocatable :: plans(:)
    type(path_support_t), allocatable :: canonical_support(:)
  contains
    procedure, public :: build => build_recoupling_system
    procedure, public :: build_plans => build_recoupling_plans
    procedure, public :: build_swaps => build_recoupling_swaps
  end type recoupling_system_t

  public :: apply_sparse_add, apply_sparse_pair_add, apply_sparse_scalar_add

contains

  subroutine build_recoupling_system(this, table, catalog)
    class(recoupling_system_t), intent(inout) :: this
    type(wigner_table_t), intent(in) :: table
    type(path_catalog_t), intent(inout) :: catalog

    call this%build_swaps(table, catalog)
    call this%build_plans(catalog)
  end subroutine build_recoupling_system

  subroutine build_recoupling_swaps(this, table, catalog)
    class(recoupling_system_t), intent(inout) :: this
    type(wigner_table_t), intent(in) :: table
    type(path_catalog_t), intent(inout) :: catalog

    integer :: length

    if (catalog%max_length < 1 .or. .not. allocated(catalog%spaces)) &
         call fail('path catalog must be built before recoupling swaps')
    this%max_length = catalog%max_length
    this%swaps_built = .false.
    this%plans_built = .false.
    if (allocated(this%swaps)) deallocate(this%swaps)
    if (allocated(this%plans)) deallocate(this%plans)
    if (allocated(this%canonical_support)) deallocate(this%canonical_support)
    allocate(this%swaps(this%max_length))

    do length = 1, this%max_length
      call build_swap_set(this%swaps(length), table, catalog%spaces(length))
      call determine_reflection_phases(catalog%spaces(length))
    end do
    this%swaps_built = .true.
  end subroutine build_recoupling_swaps

  subroutine build_recoupling_plans(this, catalog)
    class(recoupling_system_t), intent(inout) :: this
    type(path_catalog_t), intent(inout) :: catalog

    integer :: length, path

    if (.not. this%swaps_built .or. .not. allocated(this%swaps)) &
         call fail('recoupling swaps must be built before ordering plans')
    if (catalog%max_length /= this%max_length .or. &
        .not. allocated(catalog%spaces)) &
         call fail('path catalog does not match recoupling swaps')
    if (this%plans_built) return
    if (allocated(this%plans)) deallocate(this%plans)
    if (allocated(this%canonical_support)) deallocate(this%canonical_support)
    allocate(this%plans(this%max_length), &
             this%canonical_support(this%max_length))
    allocate(this%canonical_support(1)%active_path( &
         catalog%spaces(1)%number_of_paths))
    this%canonical_support(1)%active_path = &
         [(path, path=1, catalog%spaces(1)%number_of_paths)]
    do length = 2, this%max_length
      call build_ordering_plan(this%plans(length), catalog, length)
      call build_active_paths(this%plans(length), this%swaps(length), &
                              this%canonical_support)
      associate(canonical => this%plans(length)%nodes( &
                this%plans(length)%canonical_node)%active_path)
        allocate(this%canonical_support(length)%active_path(size(canonical)))
        this%canonical_support(length)%active_path = canonical
      end associate
    end do
    this%plans_built = .true.
  end subroutine build_recoupling_plans

  subroutine build_swap_set(set, table, space)
    type(swap_set_t), intent(inout) :: set
    type(wigner_table_t), intent(in) :: table
    type(path_space_t), intent(in) :: space

    integer :: position

    allocate(set%matrices(space%length))
    do position = 1, space%length
      call build_one_swap(set%matrices(position), table, space, position)
    end do
  end subroutine build_swap_set

  subroutine build_one_swap(matrix, table, space, position)
    type(sparse_matrix_t), intent(inout) :: matrix
    type(wigner_table_t), intent(in) :: table
    type(path_space_t), intent(in) :: space
    integer, intent(in) :: position

    integer, allocatable :: output_p(:), output_q(:), output_mult(:)
    integer :: block_index, entry, input_index, input_local, number_of_values
    integer :: output_index, output_local

    if (position < 1 .or. position > space%length) &
         call fail('adjacent-swap position outside closed path')
    number_of_values = 0
    do input_index = 1, space%number_of_paths
      call locate_input_block(input_index, block_index, input_local)
      do output_local = 1, table%blocks(block_index)%size
        if (abs(table%blocks(block_index)%matrix(output_local, input_local)) > &
            sparse_zero) number_of_values = number_of_values + 1
      end do
    end do

    matrix%size = space%number_of_paths
    matrix%number_of_values = number_of_values
    allocate(matrix%column_start(matrix%size+1))
    allocate(matrix%row(number_of_values), matrix%value(number_of_values))
    allocate(output_p(0:space%length), output_q(0:space%length))
    allocate(output_mult(space%length))

    entry = 1
    do input_index = 1, space%number_of_paths
      matrix%column_start(input_index) = entry
      call locate_input_block(input_index, block_index, input_local)
      do output_local = 1, table%blocks(block_index)%size
        if (abs(table%blocks(block_index)%matrix(output_local, input_local)) <= &
            sparse_zero) cycle
        output_p = space%rep_p(:, input_index)
        output_q = space%rep_q(:, input_index)
        output_mult = space%multiplicity(:, input_index)
        output_p(position) = &
             table%blocks(block_index)%paths(output_local)%middle_p
        output_q(position) = &
             table%blocks(block_index)%paths(output_local)%middle_q
        output_mult(position) = &
             table%blocks(block_index)%paths(output_local)%left_multiplicity
        if (position < space%length) then
          output_mult(position+1) = &
               table%blocks(block_index)%paths(output_local)%right_multiplicity
        else if (table%blocks(block_index)%paths(output_local)%right_multiplicity /= 0) then
          call fail('swap changed the unique adjoint-to-singlet closure')
        end if
        output_index = space%find(output_p, output_q, output_mult)
        if (output_index == 0) &
             call fail('Wigner swap produced a path absent from the catalog')
        matrix%row(entry) = output_index
        matrix%value(entry) = &
             table%blocks(block_index)%matrix(output_local, input_local)
        entry = entry + 1
      end do
    end do
    matrix%column_start(matrix%size+1) = entry
    if (entry /= number_of_values+1) call fail('sparse swap count mismatch')
    matrix%is_signed_permutation = &
         all(matrix%column_start(2:matrix%size+1)- &
             matrix%column_start(1:matrix%size) == 1) .and. &
         all(abs(abs(matrix%value)-1.0_dp) <= sparse_zero)
    if (matrix%is_signed_permutation) then
      allocate(matrix%signed_permutation(matrix%size))
      do input_index = 1, matrix%size
        entry = matrix%column_start(input_index)
        if (matrix%value(entry) > 0.0_dp) then
          matrix%signed_permutation(input_index) = matrix%row(entry)
        else
          matrix%signed_permutation(input_index) = -matrix%row(entry)
        end if
      end do
    end if

  contains

    subroutine locate_input_block(path_index, found_block, found_path)
      integer, intent(in) :: path_index
      integer, intent(out) :: found_block, found_path
      integer :: left_p, left_q, middle_p, middle_q, right_p, right_q
      integer :: left_mult, right_mult

      left_p = space%rep_p(position-1, path_index)
      left_q = space%rep_q(position-1, path_index)
      middle_p = space%rep_p(position, path_index)
      middle_q = space%rep_q(position, path_index)
      if (position < space%length) then
        right_p = space%rep_p(position+1, path_index)
        right_q = space%rep_q(position+1, path_index)
        right_mult = space%multiplicity(position+1, path_index)
      else
        right_p = 0
        right_q = 0
        right_mult = 0
      end if
      left_mult = space%multiplicity(position, path_index)
      found_block = table%find_block(left_p, left_q, right_p, right_q)
      if (found_block == 0) call fail('required outer-irrep Wigner block is absent')
      found_path = table%find_local_path(found_block, middle_p, middle_q, &
                                         left_mult, right_mult)
      if (found_path == 0) call fail('required local Wigner path is absent')
    end subroutine locate_input_block

  end subroutine build_one_swap

  subroutine determine_reflection_phases(space)
    type(path_space_t), intent(inout) :: space

    integer, allocatable :: reflected_path(:), output_p(:), output_q(:)
    integer, allocatable :: output_mult(:)
    integer :: input_path, output_path, phase_parity, position, length

    length = space%length
    allocate(reflected_path(space%number_of_paths))
    allocate(output_p(0:length), output_q(0:length), output_mult(length))
    do input_path = 1, space%number_of_paths
      ! Reversing a path conjugates and reverses its intermediate irreps and
      ! reverses its outer-multiplicity labels.  Keep constructing the image
      ! explicitly so a table with incompatible conjugation conventions is
      ! rejected even though the phase itself no longer needs any swaps.
      output_p = -1
      output_q = -1
      output_mult = -1
      output_p(0) = 0
      output_q(0) = 0
      output_p(1) = 1
      output_q(1) = 1
      output_mult(1) = 0
      do position = 2, length
        call conjugate_label(space%rep_p(length+1-position, input_path), &
                             space%rep_q(length+1-position, input_path), &
                             output_p(position), output_q(position))
        output_mult(position) = &
             space%multiplicity(length+2-position, input_path)
      end do
      output_path = space%find(output_p, output_q, output_mult)
      if (output_path == 0) call fail('reflected path is absent from catalog')
      reflected_path(input_path) = output_path

      ! In the phase convention used to generate the Wigner table, reflection
      ! of a non-self vertex has sign +1.  A self vertex R x 8 -> R has sign
      ! (-1)^multiplicity for R=8 and the opposite sign for every other R.
      ! Multiplying these local factors gives the complete path phase.
      phase_parity = 0
      do position = 1, length
        if (space%multiplicity(position, input_path) < 0) &
             call fail('negative path multiplicity in reflection phase')
        phase_parity = phase_parity + &
             space%multiplicity(position, input_path)
        if (space%rep_p(position-1, input_path) == &
            space%rep_p(position, input_path) .and. &
            space%rep_q(position-1, input_path) == &
            space%rep_q(position, input_path) .and. &
            (space%rep_p(position, input_path) /= 1 .or. &
             space%rep_q(position, input_path) /= 1)) &
             phase_parity = phase_parity + 1
      end do
      if (modulo(phase_parity, 2) == 0) then
        space%reflection_phase(input_path) = 1.0_dp
      else
        space%reflection_phase(input_path) = -1.0_dp
      end if
    end do

    do input_path = 1, space%number_of_paths
      output_path = reflected_path(input_path)
      if (reflected_path(output_path) /= input_path) &
           call fail('path reflection is not an involution')
      if (space%reflection_phase(output_path)* &
          space%reflection_phase(input_path) < 0.0_dp) &
           call fail('reflected paths have inconsistent phases')
    end do
  end subroutine determine_reflection_phases

  subroutine apply_sparse_add(matrix, source, destination, active_path, &
                              initialized)
    type(sparse_matrix_t), intent(in) :: matrix
    complex(dp), intent(in) :: source(4, matrix%size)
    complex(dp), intent(inout) :: destination(4, matrix%size)
    integer, intent(in) :: active_path(:)
    logical(c_bool), intent(inout) :: initialized(matrix%size)

    real(dp) :: coefficient
    integer :: active_index, component, entry, output_path, signed_output
    integer :: source_path

    if (matrix%is_signed_permutation) then
      do active_index = 1, size(active_path)
        source_path = active_path(active_index)
        signed_output = matrix%signed_permutation(source_path)
        output_path = abs(signed_output)
        if (.not. initialized(output_path)) then
          if (signed_output > 0) then
            destination(:, output_path) = source(:, source_path)
          else
            destination(:, output_path) = -source(:, source_path)
          end if
          initialized(output_path) = .true.
        else
          if (signed_output > 0) then
            do component = 1, 4
              destination(component, output_path) = &
                   destination(component, output_path) + &
                   source(component, source_path)
            end do
          else
            do component = 1, 4
              destination(component, output_path) = &
                   destination(component, output_path) - &
                   source(component, source_path)
            end do
          end if
        end if
      end do
      return
    end if
    do active_index = 1, size(active_path)
      source_path = active_path(active_index)
      do entry = matrix%column_start(source_path), &
                 matrix%column_start(source_path+1)-1
        output_path = matrix%row(entry)
        coefficient = matrix%value(entry)
        if (.not. initialized(output_path)) then
          do component = 1, 4
            destination(component, output_path)%re = &
                 coefficient*source(component, source_path)%re
            destination(component, output_path)%im = &
                 coefficient*source(component, source_path)%im
          end do
          initialized(output_path) = .true.
        else
          do component = 1, 4
            destination(component, output_path)%re = &
                 destination(component, output_path)%re + &
                 coefficient*source(component, source_path)%re
            destination(component, output_path)%im = &
                 destination(component, output_path)%im + &
                 coefficient*source(component, source_path)%im
          end do
        end if
      end do
    end do
  end subroutine apply_sparse_add

  subroutine apply_sparse_pair_add(matrix, source_gluon, source_tensor, &
                                   destination_gluon, destination_tensor, &
                                   active_path, initialized)
    type(sparse_matrix_t), intent(in) :: matrix
    complex(dp), intent(in) :: source_gluon(4, matrix%size)
    complex(dp), intent(in) :: source_tensor(6, matrix%size)
    complex(dp), intent(inout) :: destination_gluon(4, matrix%size)
    complex(dp), intent(inout) :: destination_tensor(6, matrix%size)
    integer, intent(in) :: active_path(:)
    logical(c_bool), intent(inout) :: initialized(matrix%size)

    real(dp) :: coefficient
    integer :: active_index, component, entry, output_path, signed_output
    integer :: source_path

    ! Lower-level vector and auxiliary-tensor currents follow the same
    ! recoupling edge.  Traverse its sparse metadata once for both currents.
    if (matrix%is_signed_permutation) then
      do active_index = 1, size(active_path)
        source_path = active_path(active_index)
        signed_output = matrix%signed_permutation(source_path)
        output_path = abs(signed_output)
        if (.not. initialized(output_path)) then
          if (signed_output > 0) then
            destination_gluon(:, output_path) = source_gluon(:, source_path)
            destination_tensor(:, output_path) = source_tensor(:, source_path)
          else
            destination_gluon(:, output_path) = -source_gluon(:, source_path)
            destination_tensor(:, output_path) = -source_tensor(:, source_path)
          end if
          initialized(output_path) = .true.
        else
          if (signed_output > 0) then
            do component = 1, 4
              destination_gluon(component, output_path) = &
                   destination_gluon(component, output_path) + &
                   source_gluon(component, source_path)
            end do
            do component = 1, 6
              destination_tensor(component, output_path) = &
                   destination_tensor(component, output_path) + &
                   source_tensor(component, source_path)
            end do
          else
            do component = 1, 4
              destination_gluon(component, output_path) = &
                   destination_gluon(component, output_path) - &
                   source_gluon(component, source_path)
            end do
            do component = 1, 6
              destination_tensor(component, output_path) = &
                   destination_tensor(component, output_path) - &
                   source_tensor(component, source_path)
            end do
          end if
        end if
      end do
      return
    end if
    do active_index = 1, size(active_path)
      source_path = active_path(active_index)
      do entry = matrix%column_start(source_path), &
                 matrix%column_start(source_path+1)-1
        output_path = matrix%row(entry)
        coefficient = matrix%value(entry)
        if (.not. initialized(output_path)) then
          do component = 1, 4
            destination_gluon(component, output_path)%re = &
                 coefficient*source_gluon(component, source_path)%re
            destination_gluon(component, output_path)%im = &
                 coefficient*source_gluon(component, source_path)%im
          end do
          do component = 1, 6
            destination_tensor(component, output_path)%re = &
                 coefficient*source_tensor(component, source_path)%re
            destination_tensor(component, output_path)%im = &
                 coefficient*source_tensor(component, source_path)%im
          end do
          initialized(output_path) = .true.
        else
          do component = 1, 4
            destination_gluon(component, output_path)%re = &
                 destination_gluon(component, output_path)%re + &
                 coefficient*source_gluon(component, source_path)%re
            destination_gluon(component, output_path)%im = &
                 destination_gluon(component, output_path)%im + &
                 coefficient*source_gluon(component, source_path)%im
          end do
          do component = 1, 6
            destination_tensor(component, output_path)%re = &
                 destination_tensor(component, output_path)%re + &
                 coefficient*source_tensor(component, source_path)%re
            destination_tensor(component, output_path)%im = &
                 destination_tensor(component, output_path)%im + &
                 coefficient*source_tensor(component, source_path)%im
          end do
        end if
      end do
    end do
  end subroutine apply_sparse_pair_add

  subroutine apply_sparse_scalar_add(matrix, source, destination, &
                                     active_path, initialized, &
                                     destination_preinitialized)
    type(sparse_matrix_t), intent(in) :: matrix
    complex(dp), intent(in) :: source(matrix%size)
    complex(dp), intent(inout) :: destination(matrix%size)
    integer, intent(in) :: active_path(:)
    logical(c_bool), intent(inout) :: initialized(matrix%size)
    logical, intent(in) :: destination_preinitialized

    real(dp) :: coefficient
    integer :: active_index, entry, output_path, signed_output, source_path

    if (destination_preinitialized) then
      if (size(active_path) == matrix%size) then
        if (matrix%is_signed_permutation) then
          do source_path = 1, matrix%size
            signed_output = matrix%signed_permutation(source_path)
            output_path = abs(signed_output)
            if (signed_output > 0) then
              destination(output_path) = destination(output_path) + &
                   source(source_path)
            else
              destination(output_path) = destination(output_path) - &
                   source(source_path)
            end if
          end do
          return
        end if
        do source_path = 1, matrix%size
          do entry = matrix%column_start(source_path), &
                     matrix%column_start(source_path+1)-1
            output_path = matrix%row(entry)
            coefficient = matrix%value(entry)
            destination(output_path)%re = destination(output_path)%re + &
                 coefficient*source(source_path)%re
            destination(output_path)%im = destination(output_path)%im + &
                 coefficient*source(source_path)%im
          end do
        end do
        return
      end if
      if (matrix%is_signed_permutation) then
        do active_index = 1, size(active_path)
          source_path = active_path(active_index)
          signed_output = matrix%signed_permutation(source_path)
          output_path = abs(signed_output)
          if (signed_output > 0) then
            destination(output_path) = destination(output_path) + &
                 source(source_path)
          else
            destination(output_path) = destination(output_path) - &
                 source(source_path)
          end if
        end do
        return
      end if
      do active_index = 1, size(active_path)
        source_path = active_path(active_index)
        do entry = matrix%column_start(source_path), &
                   matrix%column_start(source_path+1)-1
          output_path = matrix%row(entry)
          coefficient = matrix%value(entry)
          destination(output_path)%re = destination(output_path)%re + &
               coefficient*source(source_path)%re
          destination(output_path)%im = destination(output_path)%im + &
               coefficient*source(source_path)%im
        end do
      end do
      return
    end if

    if (matrix%is_signed_permutation) then
      do active_index = 1, size(active_path)
        source_path = active_path(active_index)
        signed_output = matrix%signed_permutation(source_path)
        output_path = abs(signed_output)
        if (.not. initialized(output_path)) then
          if (signed_output > 0) then
            destination(output_path) = source(source_path)
          else
            destination(output_path) = -source(source_path)
          end if
          initialized(output_path) = .true.
        else
          if (signed_output > 0) then
            destination(output_path) = destination(output_path) + &
                 source(source_path)
          else
            destination(output_path) = destination(output_path) - &
                 source(source_path)
          end if
        end if
      end do
      return
    end if
    do active_index = 1, size(active_path)
      source_path = active_path(active_index)
      do entry = matrix%column_start(source_path), &
                 matrix%column_start(source_path+1)-1
        output_path = matrix%row(entry)
        coefficient = matrix%value(entry)
        if (.not. initialized(output_path)) then
          destination(output_path)%re = coefficient*source(source_path)%re
          destination(output_path)%im = coefficient*source(source_path)%im
          initialized(output_path) = .true.
        else
          destination(output_path)%re = destination(output_path)%re + &
               coefficient*source(source_path)%re
          destination(output_path)%im = destination(output_path)%im + &
               coefficient*source(source_path)%im
        end if
      end do
    end do
  end subroutine apply_sparse_scalar_add

  subroutine build_active_paths(plan, swaps, canonical_support)
    type(recoupling_plan_t), intent(inout) :: plan
    type(swap_set_t), intent(in) :: swaps
    type(path_support_t), intent(in) :: canonical_support(:)

    logical, allocatable :: active(:, :)
    integer :: active_count, entry, index_partition, inversion, layer_entry
    integer :: left_active, left_path, node, path, position, right_active
    integer :: right_path, source_path, target_node

    allocate(active(plan%number_of_paths, plan%number_of_nodes))
    active = .false.

    do index_partition = 1, plan%number_of_partitions
      node = plan%partitions(index_partition)%initial_node
      associate(seed => plan%seeds(plan%partitions(index_partition)%left_size))
        do right_active = 1, &
             size(canonical_support(seed%right_size)%active_path)
          right_path = canonical_support(seed%right_size)% &
               active_path(right_active)
          do left_active = 1, &
               size(canonical_support(seed%left_size)%active_path)
            left_path = canonical_support(seed%left_size)% &
                 active_path(left_active)
            active(abs(seed%signed_parent_path(left_path, right_path)), node) = &
                 .true.
          end do
        end do
      end associate
    end do

    ! Sorting edges always lower the inversion count by one.  Propagate the
    ! exact structural support in that topological order; it is independent
    ! of momenta and helicities and can be reused for every event.
    do inversion = plan%max_inversions, 1, -1
      do layer_entry = plan%layer_node_start(inversion), &
                       plan%layer_node_start(inversion+1)-1
        node = plan%layer_nodes(layer_entry)
        target_node = plan%nodes(node)%next_node
        position = plan%nodes(node)%swap_position
        if (position < 1 .or. target_node < 1) &
             call fail('incomplete edge while building active paths')
        associate(matrix => swaps%matrices(position))
          plan%nodes(node)%destination_preinitialized = .true.
          do source_path = 1, matrix%size
            if (.not. active(source_path, node)) cycle
            do entry = matrix%column_start(source_path), &
                       matrix%column_start(source_path+1)-1
              if (.not. active(matrix%row(entry), target_node)) &
                   plan%nodes(node)%destination_preinitialized = .false.
              active(matrix%row(entry), target_node) = .true.
            end do
          end do
        end associate
      end do
    end do

    do node = 1, plan%number_of_nodes
      active_count = count(active(:, node))
      allocate(plan%nodes(node)%active_path(active_count))
      active_count = 0
      do path = 1, plan%number_of_paths
        if (.not. active(path, node)) cycle
        active_count = active_count+1
        plan%nodes(node)%active_path(active_count) = path
      end do
    end do
  end subroutine build_active_paths

  subroutine build_ordering_plan(plan, catalog, length)
    type(recoupling_plan_t), intent(inout) :: plan
    type(path_catalog_t), intent(in) :: catalog
    integer, intent(in) :: length

    type(ordering_node_t), allocatable :: work_nodes(:)
    integer, allocatable :: initial_order(:), layer_next(:), layer_width(:)
    integer, allocatable :: next_order(:)
    integer(i8) :: all_pattern, pattern
    integer :: current_node, initial_node, inversion, left_size, max_nodes, next_node
    integer :: node_count, partition_count, position, rank

    if (length < 2) call fail('ordering plan needs at least two external legs')
    if (length >= bit_size(pattern)-1) call fail('ordering plan exceeds bit-mask capacity')
    plan%length = length
    plan%number_of_paths = catalog%spaces(length)%number_of_paths
    all_pattern = shiftl(1_i8, length) - 1_i8
    partition_count = int(shiftl(1_i8, length-1) - 1_i8)
    plan%number_of_partitions = partition_count
    allocate(plan%partitions(partition_count), plan%seeds(length-1))
    allocate(initial_order(length+1), next_order(length+1))
    ! A partition can need at most length*(length+1)/2 adjacent swaps.
    ! Most routes coalesce much earlier, but the conservative bound keeps
    ! deeper, externally generated Wigner tables safe as well.
    max_nodes = (length*(length+1)/2+1)*partition_count + 1
    allocate(work_nodes(max_nodes))
    node_count = 0
    partition_count = 0

    do pattern = 1_i8, all_pattern-1_i8
      if (.not. btest(pattern, 0)) cycle
      partition_count = partition_count + 1
      left_size = popcnt(pattern)
      position = 0
      do rank = 1, length
        if (btest(pattern, rank-1)) then
          position = position + 1
          initial_order(position) = rank
        end if
      end do
      position = position + 1
      initial_order(position) = length+1
      do rank = length, 1, -1
        if (.not. btest(pattern, rank-1)) then
          position = position + 1
          initial_order(position) = rank
        end if
      end do
      if (position /= length+1) call fail('partition ordering has wrong length')
      initial_node = find_or_add_node(initial_order)
      plan%partitions(partition_count)%left_pattern = pattern
      plan%partitions(partition_count)%left_size = left_size
      plan%partitions(partition_count)%right_size = length-left_size
      plan%partitions(partition_count)%initial_node = initial_node

      current_node = initial_node
      do while (work_nodes(current_node)%inversions > 0)
        if (work_nodes(current_node)%next_node == 0) then
          position = first_ascending_target_inversion( &
               work_nodes(current_node)%order)
          if (position == 0) call fail('noncanonical order has no inversion')
          next_order = work_nodes(current_node)%order
          call exchange(next_order(position), next_order(position+1))
          next_node = find_or_add_node(next_order)
          work_nodes(current_node)%swap_position = position
          work_nodes(current_node)%next_node = next_node
        end if
        current_node = work_nodes(current_node)%next_node
      end do

    end do
    if (partition_count /= plan%number_of_partitions) &
         call fail('unordered partition count mismatch')

    plan%number_of_nodes = node_count
    allocate(plan%nodes(node_count))
    plan%nodes = work_nodes(1:node_count)
    plan%canonical_node = 0
    plan%max_inversions = 0
    do current_node = 1, node_count
      plan%max_inversions = max(plan%max_inversions, &
                                plan%nodes(current_node)%inversions)
      if (plan%nodes(current_node)%inversions == 0) then
        if (plan%canonical_node /= 0) call fail('multiple canonical ordering nodes')
        plan%canonical_node = current_node
      end if
    end do
    if (plan%canonical_node == 0) call fail('ordering plan lacks canonical node')

    ! Every sorting edge removes exactly one inversion.  Record compact,
    ! per-layer node and partition lists so the evaluator can reuse two
    ! adjacent-layer buffers instead of retaining one current for every DAG
    ! node.
    allocate(layer_width(0:plan%max_inversions))
    layer_width = 0
    do current_node = 1, plan%number_of_nodes
      inversion = plan%nodes(current_node)%inversions
      layer_width(inversion) = layer_width(inversion)+1
      plan%nodes(current_node)%layer_slot = layer_width(inversion)
    end do
    plan%max_layer_width = maxval(layer_width)
    allocate(plan%layer_node_start(0:plan%max_inversions+1))
    allocate(plan%layer_nodes(plan%number_of_nodes))
    plan%layer_node_start(0) = 1
    do inversion = 0, plan%max_inversions
      plan%layer_node_start(inversion+1) = &
           plan%layer_node_start(inversion)+layer_width(inversion)
    end do
    allocate(layer_next(0:plan%max_inversions))
    layer_next = plan%layer_node_start(0:plan%max_inversions)
    do current_node = 1, plan%number_of_nodes
      inversion = plan%nodes(current_node)%inversions
      plan%layer_nodes(layer_next(inversion)) = current_node
      layer_next(inversion) = layer_next(inversion)+1
    end do

    layer_width = 0
    do partition_count = 1, plan%number_of_partitions
      initial_node = plan%partitions(partition_count)%initial_node
      inversion = plan%nodes(initial_node)%inversions
      layer_width(inversion) = layer_width(inversion)+1
    end do
    allocate(plan%layer_partition_start(0:plan%max_inversions+1))
    allocate(plan%layer_partitions(plan%number_of_partitions))
    plan%layer_partition_start(0) = 1
    do inversion = 0, plan%max_inversions
      plan%layer_partition_start(inversion+1) = &
           plan%layer_partition_start(inversion)+layer_width(inversion)
    end do
    layer_next = plan%layer_partition_start(0:plan%max_inversions)
    do partition_count = 1, plan%number_of_partitions
      initial_node = plan%partitions(partition_count)%initial_node
      inversion = plan%nodes(initial_node)%inversions
      plan%layer_partitions(layer_next(inversion)) = partition_count
      layer_next(inversion) = layer_next(inversion)+1
    end do

    do left_size = 1, length-1
      call build_seed_map(plan%seeds(left_size), catalog, length, left_size)
    end do

  contains

    integer function find_or_add_node(order) result(index_node)
      integer, intent(in) :: order(:)
      integer :: candidate

      do candidate = 1, node_count
        if (all(work_nodes(candidate)%order == order)) then
          index_node = candidate
          return
        end if
      end do
      node_count = node_count + 1
      if (node_count > size(work_nodes)) call fail('ordering-node workspace overflow')
      allocate(work_nodes(node_count)%order(length+1))
      work_nodes(node_count)%order = order
      work_nodes(node_count)%inversions = count_inversions(order)
      index_node = node_count
    end function find_or_add_node

  end subroutine build_ordering_plan

  subroutine build_seed_map(seed, catalog, total_length, left_length)
    type(seed_map_t), intent(inout) :: seed
    type(path_catalog_t), intent(in) :: catalog
    integer, intent(in) :: total_length, left_length

    integer, allocatable :: parent_p(:), parent_q(:), parent_mult(:)
    integer :: child_index, left_path, parent_index, position, right_length
    integer :: right_path

    right_length = total_length-left_length
    seed%left_size = left_length
    seed%right_size = right_length
    allocate(seed%signed_parent_path( &
         catalog%spaces(left_length)%number_of_paths, &
         catalog%spaces(right_length)%number_of_paths))
    allocate(parent_p(0:total_length), parent_q(0:total_length))
    allocate(parent_mult(total_length))

    do right_path = 1, catalog%spaces(right_length)%number_of_paths
      do left_path = 1, catalog%spaces(left_length)%number_of_paths
        parent_p = -1
        parent_q = -1
        parent_mult = -1
        parent_p(0:left_length) = &
             catalog%spaces(left_length)%rep_p(:, left_path)
        parent_q(0:left_length) = &
             catalog%spaces(left_length)%rep_q(:, left_path)
        parent_mult(1:left_length) = &
             catalog%spaces(left_length)%multiplicity(:, left_path)
        position = left_length+1
        parent_p(position) = 1
        parent_q(position) = 1
        parent_mult(position) = 1
        do child_index = right_length-1, 1, -1
          position = position+1
          call conjugate_label( &
               catalog%spaces(right_length)%rep_p(child_index, right_path), &
               catalog%spaces(right_length)%rep_q(child_index, right_path), &
               parent_p(position), parent_q(position))
        end do
        position = left_length+1
        do child_index = right_length, 2, -1
          position = position+1
          parent_mult(position) = &
               catalog%spaces(right_length)%multiplicity(child_index, right_path)
        end do
        if (position /= total_length) call fail('merge seed has wrong path length')
        parent_index = catalog%spaces(total_length)%find( &
             parent_p, parent_q, parent_mult)
        if (parent_index == 0) call fail('merge seed path is absent from catalog')
        if (catalog%spaces(right_length)%reflection_phase(right_path) < &
            0.0_dp) then
          seed%signed_parent_path(left_path, right_path) = parent_index
        else
          seed%signed_parent_path(left_path, right_path) = -parent_index
        end if
      end do
    end do
  end subroutine build_seed_map

  integer function first_ascending_target_inversion(order) result(position)
    integer, intent(in) :: order(:)
    integer :: index_order

    position = 0
    do index_order = 1, size(order)-1
      if (order(index_order) > order(index_order+1)) then
        position = index_order
        return
      end if
    end do
  end function first_ascending_target_inversion

  integer function count_inversions(order) result(number)
    integer, intent(in) :: order(:)
    integer :: left, right

    number = 0
    do left = 1, size(order)-1
      do right = left+1, size(order)
        if (order(left) > order(right)) number = number+1
      end do
    end do
  end function count_inversions

  pure subroutine exchange(left, right)
    integer, intent(inout) :: left, right
    integer :: temporary

    temporary = left
    left = right
    right = temporary
  end subroutine exchange

end module recoupling_plan
