program test_recoupling
  use ampligluon_multiplet_kinds, only: dp, fail
  use multiplet_paths, only: path_catalog_t, path_space_t
  use recoupling_plan, only: recoupling_plan_t, recoupling_system_t, &
       sparse_matrix_t
  use wigner_table, only: wigner_table_t
  implicit none

  integer, parameter :: expected_counts(8) = &
       [1, 2, 8, 32, 145, 702, 3598, 19280]
  integer, parameter :: expected_layer_widths(2:8) = &
       [1, 2, 3, 4, 7, 10, 16]
  integer, parameter :: expected_positive_phases(8) = &
       [1, 1, 6, 17, 87, 373, 1954, 10056]
  integer, parameter :: expected_negative_phases(8) = &
       [0, 1, 2, 15, 58, 329, 1644, 9224]
  integer, parameter :: expected_phase_checksums(8) = &
       [1, -1, 18, 18, 2191, 16627, 595843, 8944278]
  integer, parameter :: legacy_reflection_max_length = 6
  type(path_catalog_t) :: catalog
  type(recoupling_system_t) :: recoupling
  type(wigner_table_t) :: table
  character(len=1024) :: table_file
  integer :: length, partition, seed_index

  if (command_argument_count() /= 1) then
    call fail('usage: test_recoupling WIGNER_TABLE')
  end if
  call get_command_argument(1, table_file)

  call table%load(trim(table_file))
  call catalog%build(table, 8)
  call recoupling%build(table, catalog)

  do length = 1, 8
    call require(catalog%spaces(length)%number_of_paths == expected_counts(length), &
                 'unexpected open-path count')
    call require(maxval(abs(abs(catalog%spaces(length)%reflection_phase)-1.0_dp)) &
                 < 2.0e-9_dp, 'reflection phase is not a sign')
    call check_reflection_sentinels(catalog%spaces(length), length)
    if (length <= legacy_reflection_max_length) &
         call check_legacy_reflection(catalog%spaces(length), &
                                      recoupling%swaps(length)%matrices)
    ! Expanding a length-eight swap to a dense matrix would itself require
    ! several GiB and is unrelated to the layer-workspace regression.
    if (length <= 7) call check_swaps(recoupling%swaps(length)%matrices)
  end do

  do length = 2, 8
    associate(plan => recoupling%plans(length))
      call require(plan%number_of_partitions == 2**(length-1)-1, &
                   'unexpected unordered-partition count')
      call require(plan%canonical_node > 0, 'missing canonical ordering node')
      call require(all(plan%nodes(plan%canonical_node)%order == &
                   [(seed_index, seed_index=1, length+1)]), &
                   'canonical ordering is not ascending')
      call check_layer_metadata(plan, expected_layer_widths(length))
      if (length == 8) &
           call require(plan%number_of_nodes == 248, &
                        'unexpected length-eight sorting-node count')
      do partition = 1, plan%number_of_partitions
        seed_index = plan%partitions(partition)%left_size
        call require(all(abs(plan%seeds(seed_index)%signed_parent_path) > 0 .and. &
                         abs(plan%seeds(seed_index)%signed_parent_path) <= &
                         plan%number_of_paths), &
                     'seed map contains an invalid signed parent path')
      end do
    end associate
  end do

  write(*, '(a)') 'recoupling regression: PASS'

contains

  subroutine check_layer_metadata(plan, expected_width)
    type(recoupling_plan_t), intent(in) :: plan
    integer, intent(in) :: expected_width

    logical, allocatable :: seen_node(:), seen_partition(:)
    integer :: index_entry, index_partition, inversion, node, target

    call require(plan%max_layer_width == expected_width, &
                 'unexpected maximum sorting-layer width')
    call require(lbound(plan%layer_node_start, 1) == 0 .and. &
                 ubound(plan%layer_node_start, 1) == plan%max_inversions+1, &
                 'node-layer index has the wrong bounds')
    call require(plan%layer_node_start(0) == 1 .and. &
                 plan%layer_node_start(plan%max_inversions+1) == &
                 plan%number_of_nodes+1, &
                 'node-layer index does not cover every node')

    allocate(seen_node(plan%number_of_nodes))
    seen_node = .false.
    do inversion = 0, plan%max_inversions
      call require(plan%layer_node_start(inversion+1)- &
                   plan%layer_node_start(inversion) <= &
                   plan%max_layer_width, &
                   'sorting layer exceeds the advertised width')
      do index_entry = plan%layer_node_start(inversion), &
                       plan%layer_node_start(inversion+1)-1
        node = plan%layer_nodes(index_entry)
        call require(node >= 1 .and. node <= plan%number_of_nodes, &
                     'node-layer index contains an invalid node')
        call require(.not. seen_node(node), &
                     'node-layer index contains a duplicate node')
        seen_node(node) = .true.
        call require(plan%nodes(node)%inversions == inversion, &
                     'node is stored in the wrong inversion layer')
        call require(plan%nodes(node)%layer_slot == &
                     index_entry-plan%layer_node_start(inversion)+1, &
                     'node has the wrong reusable workspace slot')
        if (inversion > 0) then
          target = plan%nodes(node)%next_node
          call require(target >= 1 .and. target <= plan%number_of_nodes, &
                       'sorting edge has an invalid target')
          call require(plan%nodes(target)%inversions == inversion-1, &
                       'sorting edge skips an inversion layer')
        end if
      end do
    end do
    call require(all(seen_node), 'node-layer index omits a sorting node')

    call require(lbound(plan%layer_partition_start, 1) == 0 .and. &
                 ubound(plan%layer_partition_start, 1) == &
                 plan%max_inversions+1, &
                 'partition-layer index has the wrong bounds')
    call require(plan%layer_partition_start(0) == 1 .and. &
                 plan%layer_partition_start(plan%max_inversions+1) == &
                 plan%number_of_partitions+1, &
                 'partition-layer index does not cover every partition')
    allocate(seen_partition(plan%number_of_partitions))
    seen_partition = .false.
    do inversion = 0, plan%max_inversions
      do index_entry = plan%layer_partition_start(inversion), &
                       plan%layer_partition_start(inversion+1)-1
        index_partition = plan%layer_partitions(index_entry)
        call require(index_partition >= 1 .and. &
                     index_partition <= plan%number_of_partitions, &
                     'partition-layer index contains an invalid partition')
        call require(.not. seen_partition(index_partition), &
                     'partition-layer index contains a duplicate partition')
        seen_partition(index_partition) = .true.
        node = plan%partitions(index_partition)%initial_node
        call require(plan%nodes(node)%inversions == inversion, &
                     'partition seed is stored in the wrong inversion layer')
      end do
    end do
    call require(all(seen_partition), &
                 'partition-layer index omits a partition seed')
  end subroutine check_layer_metadata

  subroutine check_reflection_sentinels(space, length)
    type(path_space_t), intent(in) :: space
    integer, intent(in) :: length

    integer :: checksum, path_index

    call require(count(space%reflection_phase > 0.0_dp) == &
                 expected_positive_phases(length), &
                 'unexpected number of positive reflection phases')
    call require(count(space%reflection_phase < 0.0_dp) == &
                 expected_negative_phases(length), &
                 'unexpected number of negative reflection phases')
    checksum = 0
    do path_index = 1, space%number_of_paths
      checksum = checksum + path_index*nint(space%reflection_phase(path_index))
    end do
    call require(checksum == expected_phase_checksums(length), &
                 'unexpected ordered reflection-phase checksum')
  end subroutine check_reflection_sentinels

  subroutine check_legacy_reflection(space, matrices)
    type(path_space_t), intent(in) :: space
    type(sparse_matrix_t), intent(in) :: matrices(:)

    real(dp), allocatable :: coefficients(:, :), next_coefficients(:, :)
    integer, allocatable :: current_order(:), output_p(:), output_q(:)
    integer, allocatable :: output_mult(:)
    integer :: input_path, output_path, position, index_path, length
    real(dp) :: legacy_phase, off_diagonal

    ! This is the original quadratic-memory construction retained only as a
    ! small regression oracle.  Production reflection is linear in the path
    ! count and never allocates these dense matrices.
    length = space%length
    allocate(coefficients(space%number_of_paths, space%number_of_paths))
    allocate(next_coefficients(space%number_of_paths, space%number_of_paths))
    allocate(current_order(length+1))
    coefficients = 0.0_dp
    do index_path = 1, space%number_of_paths
      coefficients(index_path, index_path) = 1.0_dp
    end do
    current_order = [(index_path, index_path=1, length+1)]

    do
      position = first_descending_target_inversion(current_order)
      if (position == 0) exit
      call apply_sparse_dense(matrices(position), coefficients, &
                              next_coefficients)
      coefficients = next_coefficients
      call exchange(current_order(position), current_order(position+1))
    end do
    call require(all(current_order == &
                 [(length+2-index_path, index_path=1, length+1)]), &
                 'legacy reflection permutation did not terminate')

    allocate(output_p(0:length), output_q(0:length), output_mult(length))
    do input_path = 1, space%number_of_paths
      output_p = -1
      output_q = -1
      output_mult = -1
      output_p(0) = 0
      output_q(0) = 0
      output_p(1) = 1
      output_q(1) = 1
      output_mult(1) = 0
      do position = 2, length
        output_p(position) = space%rep_q(length+1-position, input_path)
        output_q(position) = space%rep_p(length+1-position, input_path)
        output_mult(position) = &
             space%multiplicity(length+2-position, input_path)
      end do
      output_path = space%find(output_p, output_q, output_mult)
      call require(output_path > 0, 'legacy reflected path is absent')
      legacy_phase = coefficients(output_path, input_path)
      coefficients(output_path, input_path) = 0.0_dp
      off_diagonal = maxval(abs(coefficients(:, input_path)))
      call require(abs(legacy_phase-space%reflection_phase(input_path)) < &
                   2.0e-9_dp, 'direct and legacy reflection phases disagree')
      call require(off_diagonal < 2.0e-9_dp, &
                   'legacy reflection is not a signed permutation')
      coefficients(output_path, input_path) = legacy_phase
    end do
  end subroutine check_legacy_reflection

  integer function first_descending_target_inversion(order) result(position)
    integer, intent(in) :: order(:)
    integer :: index_order

    position = 0
    do index_order = 1, size(order)-1
      if (order(index_order) < order(index_order+1)) then
        position = index_order
        return
      end if
    end do
  end function first_descending_target_inversion

  subroutine apply_sparse_dense(matrix, input, output)
    type(sparse_matrix_t), intent(in) :: matrix
    real(dp), intent(in) :: input(:, :)
    real(dp), intent(out) :: output(:, :)
    integer :: column, entry

    call require(size(input, 1) == matrix%size .and. &
                 size(output, 1) == matrix%size .and. &
                 size(input, 2) == size(output, 2), &
                 'incompatible legacy reflection workspace')
    output = 0.0_dp
    do column = 1, matrix%size
      do entry = matrix%column_start(column), matrix%column_start(column+1)-1
        output(matrix%row(entry), :) = output(matrix%row(entry), :) + &
             matrix%value(entry)*input(column, :)
      end do
    end do
  end subroutine apply_sparse_dense

  pure subroutine exchange(left, right)
    integer, intent(inout) :: left, right
    integer :: temporary

    temporary = left
    left = right
    right = temporary
  end subroutine exchange

  subroutine check_swaps(matrices)
    type(sparse_matrix_t), intent(in) :: matrices(:)

    real(dp), allocatable :: dense(:, :), first(:), second(:)
    integer :: column, entry, index_matrix, intermediate

    do index_matrix = 1, size(matrices)
      allocate(dense(matrices(index_matrix)%size, matrices(index_matrix)%size))
      dense = 0.0_dp
      do column = 1, matrices(index_matrix)%size
        do entry = matrices(index_matrix)%column_start(column), &
                   matrices(index_matrix)%column_start(column+1)-1
          dense(matrices(index_matrix)%row(entry), column) = &
               matrices(index_matrix)%value(entry)
        end do
      end do
      call require(maxval(abs(dense-transpose(dense))) < 3.0e-9_dp, &
                   'full-path adjacent swap is not symmetric')
      allocate(first(size(dense, 1)), second(size(dense, 1)))
      do column = 1, size(dense, 2)
        first = dense(:, column)
        second = 0.0_dp
        do intermediate = 1, matrices(index_matrix)%size
          if (abs(first(intermediate)) < tiny(1.0_dp)) cycle
          do entry = matrices(index_matrix)%column_start(intermediate), &
                     matrices(index_matrix)%column_start(intermediate+1)-1
            second(matrices(index_matrix)%row(entry)) = &
                 second(matrices(index_matrix)%row(entry)) + &
                 matrices(index_matrix)%value(entry)*first(intermediate)
          end do
        end do
        second(column) = second(column)-1.0_dp
        call require(maxval(abs(second)) < 3.0e-9_dp, &
                     'full-path adjacent swap is not an involution')
      end do
      deallocate(dense, first, second)
    end do
  end subroutine check_swaps

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) call fail(message)
  end subroutine require

end program test_recoupling
