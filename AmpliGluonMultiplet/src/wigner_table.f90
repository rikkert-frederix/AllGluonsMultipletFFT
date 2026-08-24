module wigner_table
  use ampligluon_multiplet_kinds, only: dp, fail
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  character(len=*), parameter :: table_magic = 'SU3_ADJOINT_SWAP_TABLE_V1'
  real(dp), parameter :: validation_tolerance = 2.0e-10_dp

  type, public :: local_path_t
    integer :: middle_p = -1
    integer :: middle_q = -1
    integer :: left_multiplicity = -1
    integer :: right_multiplicity = -1
    integer :: left_parity = 0
    integer :: right_parity = 0
  end type local_path_t

  type, public :: swap_block_t
    integer :: left_p = -1
    integer :: left_q = -1
    integer :: right_p = -1
    integer :: right_q = -1
    integer :: size = 0
    type(local_path_t), allocatable :: paths(:)
    real(dp), allocatable :: matrix(:, :)
  end type swap_block_t

  type, public :: transition_t
    integer :: source_p = -1
    integer :: source_q = -1
    integer :: target_p = -1
    integer :: target_q = -1
    integer :: multiplicity = -1
  end type transition_t

  type, public :: wigner_table_t
    integer :: max_prefix_gluons = -1
    integer :: number_of_paths = 0
    integer :: number_of_values = 0
    type(swap_block_t), allocatable :: blocks(:)
    type(transition_t), allocatable :: transitions(:)
  contains
    procedure :: load => load_wigner_table
    procedure :: find_block
    procedure :: find_local_path
    procedure :: validate => validate_wigner_table
  end type wigner_table_t

contains

  subroutine load_wigner_table(this, filename)
    class(wigner_table_t), intent(inout) :: this
    character(len=*), intent(in) :: filename

    character(len=1024) :: line
    integer :: block_id, input_unit, ios, number_of_blocks
    integer :: path_id, left_p, left_q, right_p, right_q, block_size
    integer :: middle_p, middle_q, left_mult, right_mult, left_parity, right_parity
    integer :: output_path, input_path
    real(dp) :: coefficient
    logical :: found_blocks, found_paths, found_values

    if (allocated(this%blocks)) deallocate(this%blocks)
    if (allocated(this%transitions)) deallocate(this%transitions)
    this%max_prefix_gluons = -1
    this%number_of_paths = 0
    this%number_of_values = 0
    number_of_blocks = -1

    open(newunit=input_unit, file=filename, status='old', action='read', &
         form='formatted', iostat=ios)
    if (ios /= 0) call fail('cannot open Wigner table: '//trim(filename))

    read(input_unit, '(a)', iostat=ios) line
    if (ios /= 0 .or. trim(line) /= table_magic) then
      close(input_unit)
      call fail('invalid Wigner table magic in '//trim(filename))
    end if

    found_blocks = .false.
    do
      read(input_unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (index(line, 'MAX_PREFIX_GLUONS ') == 1) then
        read(line(len('MAX_PREFIX_GLUONS ')+1:), *, iostat=ios) &
             this%max_prefix_gluons
        if (ios /= 0) call fail('invalid MAX_PREFIX_GLUONS row')
      else if (index(line, 'NBLOCKS ') == 1) then
        read(line(len('NBLOCKS ')+1:), *, iostat=ios) number_of_blocks
        if (ios /= 0) call fail('invalid NBLOCKS row')
      else if (index(line, 'NPATHS ') == 1) then
        read(line(len('NPATHS ')+1:), *, iostat=ios) this%number_of_paths
        if (ios /= 0) call fail('invalid NPATHS row')
      else if (index(line, 'NVALUES ') == 1) then
        read(line(len('NVALUES ')+1:), *, iostat=ios) this%number_of_values
        if (ios /= 0) call fail('invalid NVALUES row')
      else if (trim(line) == 'BEGIN_BLOCKS') then
        found_blocks = .true.
        exit
      end if
    end do
    if (.not. found_blocks .or. number_of_blocks < 1) then
      close(input_unit)
      call fail('missing block metadata in Wigner table')
    end if

    allocate(this%blocks(number_of_blocks))
    do
      read(input_unit, '(a)', iostat=ios) line
      if (ios /= 0) call fail('unexpected end of Wigner block section')
      line = adjustl(line)
      if (trim(line) == 'END_BLOCKS') exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=ios) block_id, left_p, left_q, right_p, right_q, &
           block_size
      if (ios /= 0 .or. block_id < 0 .or. block_id >= number_of_blocks .or. &
          block_size < 1) call fail('invalid Wigner block row')
      this%blocks(block_id+1)%left_p = left_p
      this%blocks(block_id+1)%left_q = left_q
      this%blocks(block_id+1)%right_p = right_p
      this%blocks(block_id+1)%right_q = right_q
      this%blocks(block_id+1)%size = block_size
      allocate(this%blocks(block_id+1)%paths(block_size))
      allocate(this%blocks(block_id+1)%matrix(block_size, block_size))
      this%blocks(block_id+1)%matrix = 0.0_dp
    end do

    found_paths = .false.
    do
      read(input_unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (trim(adjustl(line)) == 'BEGIN_PATHS') then
        found_paths = .true.
        exit
      end if
    end do
    if (.not. found_paths) call fail('missing Wigner path section')

    do
      read(input_unit, '(a)', iostat=ios) line
      if (ios /= 0) call fail('unexpected end of Wigner path section')
      line = adjustl(line)
      if (trim(line) == 'END_PATHS') exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=ios) block_id, path_id, middle_p, middle_q, &
           left_mult, right_mult, left_parity, right_parity
      if (ios /= 0 .or. block_id < 0 .or. block_id >= number_of_blocks) &
           call fail('invalid Wigner path row')
      if (path_id < 0 .or. path_id >= this%blocks(block_id+1)%size) &
           call fail('Wigner path index outside its block')
      associate(path => this%blocks(block_id+1)%paths(path_id+1))
        path%middle_p = middle_p
        path%middle_q = middle_q
        path%left_multiplicity = left_mult
        path%right_multiplicity = right_mult
        path%left_parity = left_parity
        path%right_parity = right_parity
      end associate
    end do

    found_values = .false.
    do
      read(input_unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (trim(adjustl(line)) == 'BEGIN_VALUES') then
        found_values = .true.
        exit
      end if
    end do
    if (.not. found_values) call fail('missing Wigner value section')

    do
      read(input_unit, '(a)', iostat=ios) line
      if (ios /= 0) call fail('unexpected end of Wigner value section')
      line = adjustl(line)
      if (trim(line) == 'END_VALUES') exit
      if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
      read(line, *, iostat=ios) block_id, output_path, input_path, coefficient
      if (ios /= 0 .or. block_id < 0 .or. block_id >= number_of_blocks) &
           call fail('invalid Wigner value row')
      if (.not. ieee_is_finite(coefficient)) &
           call fail('non-finite coefficient in Wigner table')
      if (output_path < 0 .or. input_path < 0 .or. &
          output_path >= this%blocks(block_id+1)%size .or. &
          input_path >= this%blocks(block_id+1)%size) &
           call fail('Wigner matrix index outside its block')
      this%blocks(block_id+1)%matrix(output_path+1, input_path+1) = coefficient
    end do
    close(input_unit)

    call this%validate()
    call collect_transitions(this)
  end subroutine load_wigner_table

  integer function find_block(this, left_p, left_q, right_p, right_q) result(index_out)
    class(wigner_table_t), intent(in) :: this
    integer, intent(in) :: left_p, left_q, right_p, right_q
    integer :: index_block

    index_out = 0
    do index_block = 1, size(this%blocks)
      if (this%blocks(index_block)%left_p == left_p .and. &
          this%blocks(index_block)%left_q == left_q .and. &
          this%blocks(index_block)%right_p == right_p .and. &
          this%blocks(index_block)%right_q == right_q) then
        index_out = index_block
        return
      end if
    end do
  end function find_block

  integer function find_local_path(this, block_index, middle_p, middle_q, &
                                   left_mult, right_mult) result(index_out)
    class(wigner_table_t), intent(in) :: this
    integer, intent(in) :: block_index, middle_p, middle_q, left_mult, right_mult
    integer :: index_path

    index_out = 0
    if (block_index < 1 .or. block_index > size(this%blocks)) return
    do index_path = 1, this%blocks(block_index)%size
      associate(path => this%blocks(block_index)%paths(index_path))
        if (path%middle_p == middle_p .and. path%middle_q == middle_q .and. &
            path%left_multiplicity == left_mult .and. &
            path%right_multiplicity == right_mult) then
          index_out = index_path
          return
        end if
      end associate
    end do
  end function find_local_path

  subroutine validate_wigner_table(this)
    class(wigner_table_t), intent(in) :: this

    real(dp), allocatable :: gram(:, :)
    real(dp) :: error
    integer :: block_index, index_path

    if (.not. allocated(this%blocks)) call fail('Wigner table is not allocated')
    if (this%max_prefix_gluons < 0) call fail('invalid Wigner table cutoff')
    do block_index = 1, size(this%blocks)
      associate(block => this%blocks(block_index))
        if (block%size < 1 .or. .not. allocated(block%paths) .or. &
            .not. allocated(block%matrix)) call fail('incomplete Wigner block')
        if (.not. all(ieee_is_finite(block%matrix))) &
             call fail('non-finite coefficient in Wigner table')
        do index_path = 1, block%size
          if (block%paths(index_path)%middle_p < 0 .or. &
              block%paths(index_path)%middle_q < 0 .or. &
              block%paths(index_path)%left_multiplicity < 0 .or. &
              block%paths(index_path)%right_multiplicity < 0) &
               call fail('incomplete local path in Wigner table')
        end do
        error = maxval(abs(block%matrix - transpose(block%matrix)))
        if (error > validation_tolerance) &
             call fail('non-symmetric Wigner swap block')
        gram = matmul(transpose(block%matrix), block%matrix)
        do index_path = 1, block%size
          gram(index_path, index_path) = gram(index_path, index_path) - 1.0_dp
        end do
        error = maxval(abs(gram))
        if (error > validation_tolerance*real(max(1, block%size), dp)) &
             call fail('non-orthogonal Wigner swap block')
      end associate
    end do
  end subroutine validate_wigner_table

  subroutine collect_transitions(this)
    class(wigner_table_t), intent(inout) :: this

    type(transition_t), allocatable :: work(:)
    integer :: block_index, count, path_index

    allocate(work(2*this%number_of_paths))
    count = 0
    do block_index = 1, size(this%blocks)
      do path_index = 1, this%blocks(block_index)%size
        associate(block => this%blocks(block_index), &
                  path => this%blocks(block_index)%paths(path_index))
          call add_unique(block%left_p, block%left_q, path%middle_p, &
                          path%middle_q, path%left_multiplicity)
          call add_unique(path%middle_p, path%middle_q, block%right_p, &
                          block%right_q, path%right_multiplicity)
        end associate
      end do
    end do
    allocate(this%transitions(count))
    this%transitions = work(1:count)

  contains

    subroutine add_unique(source_p, source_q, target_p, target_q, multiplicity)
      integer, intent(in) :: source_p, source_q, target_p, target_q, multiplicity
      integer :: index_transition

      do index_transition = 1, count
        if (work(index_transition)%source_p == source_p .and. &
            work(index_transition)%source_q == source_q .and. &
            work(index_transition)%target_p == target_p .and. &
            work(index_transition)%target_q == target_q .and. &
            work(index_transition)%multiplicity == multiplicity) return
      end do
      count = count + 1
      if (count > size(work)) call fail('transition workspace overflow')
      work(count)%source_p = source_p
      work(count)%source_q = source_q
      work(count)%target_p = target_p
      work(count)%target_q = target_q
      work(count)%multiplicity = multiplicity
    end subroutine add_unique

  end subroutine collect_transitions

end module wigner_table
