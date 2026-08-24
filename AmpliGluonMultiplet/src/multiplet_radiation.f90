module multiplet_radiation
  use ampligluon_multiplet_kinds, only: dp, fail
  use multiplet_paths, only: path_catalog_t
  use recoupling_plan, only: recoupling_system_t, sparse_matrix_t
  implicit none
  private

  type :: radiation_level_t
    integer :: source_size = 0
    integer :: target_size = 0
    integer, allocatable :: appended_path(:)
  end type radiation_level_t

  ! Work buffers are deliberately separate from the immutable radiation
  ! system.  One workspace can be retained across events and grows only when
  ! a larger path space is first encountered.  The evaluator remains serial;
  ! independent callers can use independent workspaces if required.
  type, public :: multiplet_radiation_workspace_t
    private
    complex(dp), allocatable :: source_first(:)
    complex(dp), allocatable :: source_second(:)
    complex(dp), allocatable :: target_first(:)
    complex(dp), allocatable :: target_second(:)
    logical :: result_in_first = .true.
    integer :: result_size = 0
  end type multiplet_radiation_workspace_t

  ! A level maps a closed basis with old_length open vertices (and therefore
  ! old_length+1 gluons) to the next closed multiplet basis.  For emitter i,
  ! it moves i to the old closure, appends the normalized antisymmetric
  ! 8 x 8 -> 8 vertex (outer multiplicity one), and undoes the target swaps.
  ! All intermediate and final vectors use the existing orthonormal path
  ! bases; no trace, DDM, or adjoint-index representation is introduced.
  type, public :: multiplet_radiation_system_t
    private
    integer :: max_old_length = 0
    type(radiation_level_t), allocatable :: levels(:)
  contains
    procedure :: build => build_radiation_system
    procedure :: apply => apply_radiation
    procedure :: add => add_radiation
    procedure :: number_of_emitters => radiation_number_of_emitters
    procedure :: source_dimension => radiation_source_dimension
    procedure :: target_dimension => radiation_target_dimension
  end type multiplet_radiation_system_t

contains

  subroutine build_radiation_system(this, catalog)
    class(multiplet_radiation_system_t), intent(inout) :: this
    type(path_catalog_t), intent(in) :: catalog

    logical, allocatable :: target_seen(:)
    integer, allocatable :: rep_p(:), rep_q(:), multiplicity(:)
    integer :: new_length, old_length, path, target_path

    if (.not. allocated(catalog%spaces) .or. catalog%max_length < 2) &
         call fail('radiation system needs path spaces through length two')
    if (size(catalog%spaces) < catalog%max_length) &
         call fail('incomplete path catalog in radiation system')

    if (allocated(this%levels)) deallocate(this%levels)
    this%max_old_length = catalog%max_length-1
    allocate(this%levels(this%max_old_length))

    do old_length = 1, this%max_old_length
      new_length = old_length+1
      associate(level => this%levels(old_length), &
                source => catalog%spaces(old_length), &
                target => catalog%spaces(new_length))
        if (source%length /= old_length .or. target%length /= new_length) &
             call fail('misindexed path catalog in radiation system')
        level%source_size = source%number_of_paths
        level%target_size = target%number_of_paths
        if (level%source_size < 1 .or. level%target_size < 1) &
             call fail('empty path space in radiation system')
        allocate(level%appended_path(level%source_size))
        allocate(rep_p(0:new_length), rep_q(0:new_length), &
                 multiplicity(new_length), target_seen(level%target_size))
        target_seen = .false.

        do path = 1, level%source_size
          rep_p = -1
          rep_q = -1
          multiplicity = -1
          rep_p(0:old_length) = source%rep_p(:, path)
          rep_q(0:old_length) = source%rep_q(:, path)
          multiplicity(1:old_length) = source%multiplicity(:, path)
          rep_p(new_length) = 1
          rep_q(new_length) = 1
          multiplicity(new_length) = 1
          target_path = target%find(rep_p, rep_q, multiplicity)
          if (target_path == 0) &
               call fail('antisymmetric radiation path is absent from catalog')
          if (target_seen(target_path)) &
               call fail('antisymmetric radiation append map is not injective')
          level%appended_path(path) = target_path
          target_seen(target_path) = .true.
        end do

        deallocate(rep_p, rep_q, multiplicity, target_seen)
      end associate
    end do
  end subroutine build_radiation_system

  subroutine apply_radiation(this, recoupling, old_length, emitter, source, &
                             destination, workspace)
    class(multiplet_radiation_system_t), intent(in) :: this
    type(recoupling_system_t), intent(in) :: recoupling
    integer, intent(in) :: old_length, emitter
    complex(dp), intent(in) :: source(:)
    complex(dp), intent(out) :: destination(:)
    type(multiplet_radiation_workspace_t), intent(inout) :: workspace

    call construct_radiation(this, recoupling, old_length, emitter, source, &
                             workspace)
    if (size(destination) /= workspace%result_size) &
         call fail('radiation destination has wrong size')
    if (workspace%result_in_first) then
      destination = workspace%target_first(1:workspace%result_size)
    else
      destination = workspace%target_second(1:workspace%result_size)
    end if
  end subroutine apply_radiation

  subroutine add_radiation(this, recoupling, old_length, emitter, source, &
                           weight, destination, workspace)
    class(multiplet_radiation_system_t), intent(in) :: this
    type(recoupling_system_t), intent(in) :: recoupling
    integer, intent(in) :: old_length, emitter
    complex(dp), intent(in) :: source(:)
    complex(dp), intent(in) :: weight
    complex(dp), intent(inout) :: destination(:)
    type(multiplet_radiation_workspace_t), intent(inout) :: workspace

    call construct_radiation(this, recoupling, old_length, emitter, source, &
                             workspace)
    if (size(destination) /= workspace%result_size) &
         call fail('radiation accumulation destination has wrong size')
    if (workspace%result_in_first) then
      destination = destination + &
           weight*workspace%target_first(1:workspace%result_size)
    else
      destination = destination + &
           weight*workspace%target_second(1:workspace%result_size)
    end if
  end subroutine add_radiation

  subroutine construct_radiation(this, recoupling, old_length, emitter, &
                                 source, workspace)
    class(multiplet_radiation_system_t), intent(in) :: this
    type(recoupling_system_t), intent(in) :: recoupling
    integer, intent(in) :: old_length, emitter
    complex(dp), intent(in) :: source(:)
    type(multiplet_radiation_workspace_t), intent(inout) :: workspace

    logical :: source_in_first, target_in_first
    integer :: path, position, source_size, target_size

    call validate_level(this, old_length)
    associate(level => this%levels(old_length))
      source_size = level%source_size
      target_size = level%target_size
      if (size(source) /= source_size) &
           call fail('radiation source has wrong size')
      if (emitter < 1 .or. emitter > old_length+1) &
           call fail('radiation emitter is outside the old closed basis')
      if (.not. allocated(recoupling%swaps) .or. &
          recoupling%max_length < old_length+1) &
           call fail('radiation needs adjacent swaps through its target length')
      if (.not. allocated(recoupling%swaps(old_length)%matrices) .or. &
          .not. allocated(recoupling%swaps(old_length+1)%matrices)) &
           call fail('uninitialized adjacent swaps in radiation system')
      if (size(recoupling%swaps(old_length)%matrices) /= old_length .or. &
          size(recoupling%swaps(old_length+1)%matrices) /= old_length+1) &
           call fail('incomplete adjacent swaps in radiation system')

      call reserve_workspace(workspace, source_size, target_size)
      workspace%source_first(1:source_size) = source
      source_in_first = .true.

      ! Put the emitting old gluon in the closing position.  The last old
      ! swap exchanges the last open gluon with that closure.
      do position = emitter, old_length
        if (source_in_first) then
          call apply_swap(recoupling%swaps(old_length)%matrices(position), &
                          workspace%source_first(1:source_size), &
                          workspace%source_second(1:source_size))
        else
          call apply_swap(recoupling%swaps(old_length)%matrices(position), &
                          workspace%source_second(1:source_size), &
                          workspace%source_first(1:source_size))
        end if
        source_in_first = .not. source_in_first
      end do

      workspace%target_first(1:target_size) = (0.0_dp, 0.0_dp)
      if (source_in_first) then
        do path = 1, source_size
          workspace%target_first(level%appended_path(path)) = &
               workspace%source_first(path)
        end do
      else
        do path = 1, source_size
          workspace%target_first(level%appended_path(path)) = &
               workspace%source_second(path)
        end do
      end if
      target_in_first = .true.

      ! Express every emitter contribution in the same target convention.
      ! These are the inverse of the old emitter-to-closure swaps embedded in
      ! the new path space.  Wigner swaps are symmetric involutions, hence the
      ! reversed positional order is their inverse.
      do position = old_length, emitter, -1
        if (target_in_first) then
          call apply_swap(recoupling%swaps(old_length+1)%matrices(position), &
                          workspace%target_first(1:target_size), &
                          workspace%target_second(1:target_size))
        else
          call apply_swap(recoupling%swaps(old_length+1)%matrices(position), &
                          workspace%target_second(1:target_size), &
                          workspace%target_first(1:target_size))
        end if
        target_in_first = .not. target_in_first
      end do

      workspace%result_in_first = target_in_first
      workspace%result_size = target_size
    end associate
  end subroutine construct_radiation

  subroutine reserve_workspace(workspace, source_size, target_size)
    type(multiplet_radiation_workspace_t), intent(inout) :: workspace
    integer, intent(in) :: source_size, target_size

    if (.not. allocated(workspace%source_first)) then
      allocate(workspace%source_first(source_size), &
               workspace%source_second(source_size))
    else if (size(workspace%source_first) < source_size) then
      deallocate(workspace%source_first, workspace%source_second)
      allocate(workspace%source_first(source_size), &
               workspace%source_second(source_size))
    end if

    if (.not. allocated(workspace%target_first)) then
      allocate(workspace%target_first(target_size), &
               workspace%target_second(target_size))
    else if (size(workspace%target_first) < target_size) then
      deallocate(workspace%target_first, workspace%target_second)
      allocate(workspace%target_first(target_size), &
               workspace%target_second(target_size))
    end if
  end subroutine reserve_workspace

  subroutine apply_swap(matrix, source, destination)
    type(sparse_matrix_t), intent(in) :: matrix
    complex(dp), intent(in) :: source(:)
    complex(dp), intent(out) :: destination(:)

    integer :: column, entry, signed_row

    if (size(source) /= matrix%size .or. size(destination) /= matrix%size) &
         call fail('radiation swap vector has wrong size')
    if (matrix%is_signed_permutation) then
      do column = 1, matrix%size
        signed_row = matrix%signed_permutation(column)
        if (signed_row > 0) then
          destination(signed_row) = source(column)
        else
          destination(-signed_row) = -source(column)
        end if
      end do
      return
    end if

    destination = (0.0_dp, 0.0_dp)
    do column = 1, matrix%size
      do entry = matrix%column_start(column), matrix%column_start(column+1)-1
        destination(matrix%row(entry)) = destination(matrix%row(entry)) + &
             matrix%value(entry)*source(column)
      end do
    end do
  end subroutine apply_swap

  integer function radiation_number_of_emitters(this, old_length) result(count)
    class(multiplet_radiation_system_t), intent(in) :: this
    integer, intent(in) :: old_length

    call validate_level(this, old_length)
    count = old_length+1
  end function radiation_number_of_emitters

  integer function radiation_source_dimension(this, old_length) result(count)
    class(multiplet_radiation_system_t), intent(in) :: this
    integer, intent(in) :: old_length

    call validate_level(this, old_length)
    count = this%levels(old_length)%source_size
  end function radiation_source_dimension

  integer function radiation_target_dimension(this, old_length) result(count)
    class(multiplet_radiation_system_t), intent(in) :: this
    integer, intent(in) :: old_length

    call validate_level(this, old_length)
    count = this%levels(old_length)%target_size
  end function radiation_target_dimension

  subroutine validate_level(this, old_length)
    class(multiplet_radiation_system_t), intent(in) :: this
    integer, intent(in) :: old_length

    if (.not. allocated(this%levels) .or. &
        old_length < 1 .or. old_length > this%max_old_length) &
         call fail('radiation level is not initialized')
  end subroutine validate_level

end module multiplet_radiation
