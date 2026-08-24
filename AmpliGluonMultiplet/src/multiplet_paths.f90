module multiplet_paths
  use ampligluon_multiplet_kinds, only: dp, fail, i8
  use wigner_table, only: transition_t, wigner_table_t
  implicit none
  private

  type, public :: path_space_t
    integer :: length = 0
    integer :: number_of_paths = 0
    integer, allocatable :: rep_p(:, :)
    integer, allocatable :: rep_q(:, :)
    integer, allocatable :: multiplicity(:, :)
    real(dp), allocatable :: reflection_phase(:)
    integer :: lookup_size = 0
    integer, allocatable :: lookup(:)
  contains
    procedure :: find => find_path
  end type path_space_t

  type, public :: path_catalog_t
    integer :: max_length = 0
    type(path_space_t), allocatable :: spaces(:)
  contains
    procedure :: build => build_path_catalog
  end type path_catalog_t

  public :: conjugate_label

contains

  pure subroutine conjugate_label(input_p, input_q, output_p, output_q)
    integer, intent(in) :: input_p, input_q
    integer, intent(out) :: output_p, output_q

    output_p = input_q
    output_q = input_p
  end subroutine conjugate_label

  subroutine build_path_catalog(this, table, max_length)
    class(path_catalog_t), intent(inout) :: this
    type(wigner_table_t), intent(in) :: table
    integer, intent(in) :: max_length

    integer :: length

    if (max_length < 1) call fail('path catalog needs a positive length')
    if (.not. allocated(table%transitions)) &
         call fail('cannot build paths before loading the Wigner table')
    if (allocated(this%spaces)) deallocate(this%spaces)
    this%max_length = max_length
    allocate(this%spaces(max_length))
    do length = 1, max_length
      call build_path_space(this%spaces(length), table%transitions, length)
    end do
  end subroutine build_path_catalog

  subroutine build_path_space(space, transitions, length)
    type(path_space_t), intent(inout) :: space
    type(transition_t), intent(in) :: transitions(:)
    integer, intent(in) :: length

    integer, allocatable :: work_p(:), work_q(:), work_multiplicity(:)
    integer :: count

    allocate(work_p(0:length), work_q(0:length), work_multiplicity(length))
    work_p = -1
    work_q = -1
    work_multiplicity = -1
    work_p(0) = 0
    work_q(0) = 0
    count = 0
    call walk_count(0, 0, 0)
    if (count < 1) call fail('no open adjoint paths at requested length')

    space%length = length
    space%number_of_paths = count
    allocate(space%rep_p(0:length, count), space%rep_q(0:length, count))
    allocate(space%multiplicity(length, count))
    allocate(space%reflection_phase(count))
    space%reflection_phase = 0.0_dp
    count = 0
    call walk_fill(0, 0, 0)
    if (count /= space%number_of_paths) call fail('path enumeration changed between passes')
    call build_path_lookup(space)

  contains

    recursive subroutine walk_count(step, source_p, source_q)
      integer, intent(in) :: step, source_p, source_q
      integer :: index_transition

      if (step == length) then
        if (source_p == 1 .and. source_q == 1) count = count + 1
        return
      end if
      do index_transition = 1, size(transitions)
        if (transitions(index_transition)%source_p /= source_p .or. &
            transitions(index_transition)%source_q /= source_q) cycle
        call walk_count(step+1, transitions(index_transition)%target_p, &
                       transitions(index_transition)%target_q)
      end do
    end subroutine walk_count

    recursive subroutine walk_fill(step, source_p, source_q)
      integer, intent(in) :: step, source_p, source_q
      integer :: index_transition

      if (step == length) then
        if (source_p == 1 .and. source_q == 1) then
          count = count + 1
          space%rep_p(:, count) = work_p
          space%rep_q(:, count) = work_q
          space%multiplicity(:, count) = work_multiplicity
        end if
        return
      end if
      do index_transition = 1, size(transitions)
        if (transitions(index_transition)%source_p /= source_p .or. &
            transitions(index_transition)%source_q /= source_q) cycle
        work_p(step+1) = transitions(index_transition)%target_p
        work_q(step+1) = transitions(index_transition)%target_q
        work_multiplicity(step+1) = transitions(index_transition)%multiplicity
        call walk_fill(step+1, transitions(index_transition)%target_p, &
                       transitions(index_transition)%target_q)
      end do
    end subroutine walk_fill

  end subroutine build_path_space

  subroutine build_path_lookup(space)
    type(path_space_t), intent(inout) :: space

    integer(i8) :: hash
    integer :: bucket, path_index

    space%lookup_size = 1
    do while (space%lookup_size < 2*space%number_of_paths)
      if (space%lookup_size > shiftr(huge(space%lookup_size), 1)) &
           call fail('path lookup size overflow')
      space%lookup_size = 2*space%lookup_size
    end do
    if (allocated(space%lookup)) deallocate(space%lookup)
    allocate(space%lookup(space%lookup_size))
    space%lookup = 0
    do path_index = 1, space%number_of_paths
      hash = path_hash(space%rep_p(:, path_index), &
                       space%rep_q(:, path_index), &
                       space%multiplicity(:, path_index))
      bucket = 1+int(iand(hash, int(space%lookup_size-1, i8)))
      do while (space%lookup(bucket) /= 0)
        if (same_path(space, space%lookup(bucket), &
                      space%rep_p(:, path_index), &
                      space%rep_q(:, path_index), &
                      space%multiplicity(:, path_index))) &
             call fail('duplicate path in catalog')
        bucket = modulo(bucket, space%lookup_size)+1
      end do
      space%lookup(bucket) = path_index
    end do
  end subroutine build_path_lookup

  integer function find_path(this, rep_p, rep_q, multiplicity) result(index_out)
    class(path_space_t), intent(in) :: this
    integer, intent(in) :: rep_p(0:), rep_q(0:), multiplicity(:)
    integer(i8) :: hash
    integer :: bucket, path_index, probes

    index_out = 0
    if (ubound(rep_p, 1) /= this%length .or. &
        ubound(rep_q, 1) /= this%length .or. &
        size(multiplicity) /= this%length) return
    if (.not. allocated(this%lookup) .or. this%lookup_size < 1) &
         call fail('path lookup is not initialized')
    hash = path_hash(rep_p, rep_q, multiplicity)
    bucket = 1+int(iand(hash, int(this%lookup_size-1, i8)))
    do probes = 1, this%lookup_size
      path_index = this%lookup(bucket)
      if (path_index == 0) return
      if (same_path(this, path_index, rep_p, rep_q, multiplicity)) then
        index_out = path_index
        return
      end if
      bucket = modulo(bucket, this%lookup_size)+1
    end do
  end function find_path

  pure logical function same_path(space, path_index, rep_p, rep_q, &
                                  multiplicity) result(matches)
    type(path_space_t), intent(in) :: space
    integer, intent(in) :: path_index
    integer, intent(in) :: rep_p(0:), rep_q(0:), multiplicity(:)

    matches = all(space%rep_p(:, path_index) == rep_p) .and. &
              all(space%rep_q(:, path_index) == rep_q) .and. &
              all(space%multiplicity(:, path_index) == multiplicity)
  end function same_path

  pure integer(i8) function path_hash(rep_p, rep_q, multiplicity) result(hash)
    integer, intent(in) :: rep_p(0:), rep_q(0:), multiplicity(:)
    integer :: position

    hash = 1469598103934665603_i8
    call mix_hash(hash, rep_p(0)+17)
    call mix_hash(hash, rep_q(0)+37)
    do position = 1, ubound(rep_p, 1)
      call mix_hash(hash, rep_p(position)+17)
      call mix_hash(hash, rep_q(position)+37)
      call mix_hash(hash, multiplicity(position)+53)
    end do
  end function path_hash

  pure subroutine mix_hash(hash, value)
    integer(i8), intent(inout) :: hash
    integer, intent(in) :: value

    hash = ieor(hash, int(value, i8))
    hash = ieor(hash, shiftl(hash, 13))
    hash = ieor(hash, shiftr(hash, 7))
    hash = ieor(hash, shiftl(hash, 17))
  end subroutine mix_hash

end module multiplet_paths
