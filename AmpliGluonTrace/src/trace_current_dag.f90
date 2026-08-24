module trace_current_dag
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ampligluon_common, only: dp, fail, i64
  use gluon_kinematics, only: aux_tensor_gluon_terminal, &
       aux_tensor_gluon_to_gluon, gluon_aux_tensor_terminal, &
       gluon_aux_tensor_to_gluon, massless_vector_propagator, three_gluon, &
       three_gluon_terminal, two_gluon_to_aux_tensor
  implicit none
  private

  integer, parameter :: interaction_three_gluon = 1
  integer, parameter :: interaction_tensor = 2
  integer, parameter :: interaction_tensor_gluon = 3
  integer, parameter :: interaction_gluon_tensor = 4

  type :: trace_current_level_t
    integer :: number_of_currents = 0
    integer, allocatable :: subset_mask(:)
    complex(dp), allocatable :: vector(:, :)
    complex(dp), allocatable :: tensor(:, :)
  end type trace_current_level_t

  type :: trace_interaction_level_t
    ! Every allowed vertex kind for one ordered child pair has the same signed
    ! parent targets.  Store the pair and target list once; evaluation combines
    ! the vector-valued rules before scattering their shared contribution.
    integer :: number_of_interactions = 0
    integer :: number_of_attachments = 0
    integer :: number_of_pairs = 0
    integer :: number_of_pair_attachments = 0
    integer, allocatable :: left_current(:)
    integer, allocatable :: right_current(:)
    integer, allocatable :: left_length_start(:)
    integer, allocatable :: attachment_start(:)
    integer, allocatable :: attachment_target(:)
  end type trace_interaction_level_t

  type :: trace_word_level_t
    integer, allocatable :: word(:, :)
    integer, allocatable :: maximum(:)
  end type trace_word_level_t

  type, public :: trace_current_dag_t
    private
    logical :: initialized = .false.
    logical :: fixed_first_value = .false.
    integer :: degree_value = 0
    integer :: number_of_orders_value = 0
    type(trace_current_level_t), allocatable :: level(:)
    type(trace_interaction_level_t), allocatable :: interactions(:)
    integer, allocatable :: signed_order_current(:)
    real(dp), allocatable :: subset_momentum(:, :)
  contains
    procedure, public :: evaluate => evaluate_trace_current_dag
    procedure, public :: initialize => initialize_trace_current_dag
    procedure, public :: number_of_attachments
    procedure, public :: number_of_interactions
    procedure, public :: number_of_orders
    procedure, public :: number_of_word_slots
  end type trace_current_dag_t

contains

  subroutine initialize_trace_current_dag(this, degree, fixed_first)
    class(trace_current_dag_t), intent(inout) :: this
    integer, intent(in) :: degree
    logical, intent(in), optional :: fixed_first

    integer(i64) :: number_of_orders_64, rank_map_size_64
    integer :: length, number_of_subsets
    logical :: use_fixed_first
    type(trace_word_level_t), allocatable :: word_level(:)

    if (degree < 3) call fail('trace current DAG degree must be at least three')
    if (degree >= bit_size(0)-1) &
         call fail('too many trace-current legs for subset masks')
    use_fixed_first = .false.
    if (present(fixed_first)) use_fixed_first = fixed_first
    rank_map_size_64 = falling_factorial_i64(degree, degree)
    if (rank_map_size_64 > int(huge(0), i64)) &
         call fail('trace current DAG rank map exceeds array-index capacity')
    if (use_fixed_first) then
      number_of_orders_64 = falling_factorial_i64(degree-1, degree-1)
    else
      number_of_orders_64 = rank_map_size_64
    end if
    if (number_of_orders_64 > int(huge(0), i64)) &
         call fail('trace current DAG exceeds array-index capacity')
    call ensure_schedule_fits_default_integer(degree, use_fixed_first)

    this%initialized = .false.
    this%fixed_first_value = use_fixed_first
    this%degree_value = degree
    this%number_of_orders_value = int(number_of_orders_64)
    if (allocated(this%level)) deallocate(this%level)
    if (allocated(this%interactions)) deallocate(this%interactions)
    if (allocated(this%signed_order_current)) &
         deallocate(this%signed_order_current)
    if (allocated(this%subset_momentum)) deallocate(this%subset_momentum)
    allocate(this%level(degree), this%interactions(2:degree))
    allocate(word_level(degree))

    do length = 1, degree
      call enumerate_canonical_words(degree, length, use_fixed_first, &
           word_level(length)%word)
      allocate(word_level(length)%maximum( &
           size(word_level(length)%word, 2)))
      word_level(length)%maximum = maxval(word_level(length)%word, dim=1)
      this%level(length)%number_of_currents = &
           size(word_level(length)%word, 2)
      allocate(this%level(length)%subset_mask( &
           size(word_level(length)%word, 2)))
      call fill_subset_masks(word_level(length)%word, &
           this%level(length)%subset_mask)
      if (length == degree) then
        ! Terminal currents are contracted with the final wavefunction while
        ! they are accumulated, so only one complex scalar is stored here.
        allocate(this%level(length)%vector(1, &
             size(word_level(length)%word, 2)))
      else
        allocate(this%level(length)%vector(4, &
             size(word_level(length)%word, 2)))
      end if
      if (length >= 2 .and. length < degree) &
           allocate(this%level(length)%tensor(6, &
                size(word_level(length)%word, 2)))
    end do

    do length = 2, degree
      call build_interaction_level(this, length, word_level)
    end do
    call build_order_map(this, word_level(degree)%word)
    number_of_subsets = shiftl(1, degree)
    allocate(this%subset_momentum(0:3, 0:number_of_subsets-1))
    this%initialized = .true.
  end subroutine initialize_trace_current_dag

  subroutine build_interaction_level(this, parent_length, word_level)
    class(trace_current_dag_t), intent(inout) :: this
    integer, intent(in) :: parent_length
    type(trace_word_level_t), intent(in) :: word_level(:)

    integer :: first, first_length, kinds, pair
    integer :: number, number_of_attachments, number_of_pairs
    integer :: number_of_pair_attachments
    integer :: cursor, second, second_length, targets
    integer, allocatable :: parent_map(:)
    integer(i64), allocatable :: rank_weight(:)

    call count_interaction_level(this, parent_length, word_level, number, &
         number_of_attachments, number_of_pairs, &
         number_of_pair_attachments)
    this%interactions(parent_length)%number_of_interactions = number
    this%interactions(parent_length)%number_of_attachments = &
         number_of_attachments
    this%interactions(parent_length)%number_of_pairs = number_of_pairs
    this%interactions(parent_length)%number_of_pair_attachments = &
         number_of_pair_attachments
    allocate(this%interactions(parent_length)%left_current(number_of_pairs))
    allocate(this%interactions(parent_length)%right_current(number_of_pairs))
    allocate(this%interactions(parent_length)% &
         left_length_start(0:parent_length-1))
    this%interactions(parent_length)%left_length_start = 0
    allocate(this%interactions(parent_length)% &
         attachment_start(0:number_of_pairs))
    this%interactions(parent_length)%attachment_start = 0

    pair = 0
    targets = 0
    do first_length = 1, parent_length-1
      second_length = parent_length-first_length
      if (.not. this%fixed_first_value) &
           targets = full_interaction_target_count(first_length, second_length)
      do first = 1, size(word_level(first_length)%word, 2)
        do second = 1, size(word_level(second_length)%word, 2)
          if (iand(this%level(first_length)%subset_mask(first), &
               this%level(second_length)%subset_mask(second)) /= 0) cycle
          if (word_level(first_length)%maximum(first) >= &
               word_level(second_length)%maximum(second)) cycle
          if (this%fixed_first_value) targets = count_interaction_targets( &
               word_level(first_length)%word(:, first), &
               word_level(second_length)%word(:, second), .true.)
          pair = pair+1
          call set_interaction_pair(this%interactions(parent_length), &
               pair, first, second)
          this%interactions(parent_length)%attachment_start(pair) = targets
        end do
      end do
      this%interactions(parent_length)%left_length_start(first_length) = pair
    end do
    if (pair /= number_of_pairs) &
         call fail('trace current DAG interaction-pair count mismatch')
    if (sum(this%interactions(parent_length)% &
         attachment_start(1:number_of_pairs)) /= &
         number_of_pair_attachments) &
         call fail('trace current DAG pair-attachment count mismatch')

    do pair = 1, number_of_pairs
      this%interactions(parent_length)%attachment_start(pair) = &
           this%interactions(parent_length)%attachment_start(pair-1)+ &
           this%interactions(parent_length)%attachment_start(pair)
    end do
    allocate(this%interactions(parent_length)%attachment_target( &
         number_of_pair_attachments))
    allocate(rank_weight(parent_length))
    call fill_partial_permutation_weights(this%degree_value, parent_length, &
         rank_weight)
    call build_signed_word_map(this%degree_value, parent_length, &
         this%fixed_first_value, word_level(parent_length)%word, parent_map)
    do first_length = 1, parent_length-1
      second_length = parent_length-first_length
      do pair = this%interactions(parent_length)% &
           left_length_start(first_length-1)+1, &
           this%interactions(parent_length)%left_length_start(first_length)
        first = this%interactions(parent_length)%left_current(pair)
        second = this%interactions(parent_length)%right_current(pair)
        cursor = this%interactions(parent_length)% &
             attachment_start(pair-1)
        call append_interaction_targets( &
             this%interactions(parent_length), &
             word_level(first_length)%word(:, first), &
             word_level(second_length)%word(:, second), parent_map, &
             this%degree_value, rank_weight, cursor)
        if (cursor /= this%interactions(parent_length)% &
             attachment_start(pair)) &
             call fail('trace current DAG attachment cursor mismatch')
      end do
    end do

    ! Retain the public expanded schedule counts even though each child pair
    ! and its signed target list are stored and traversed only once.
    number = 0
    number_of_attachments = 0
    do first_length = 1, parent_length-1
      second_length = parent_length-first_length
      kinds = number_of_interaction_kinds(first_length, second_length, &
           parent_length, this%degree_value)
      do pair = this%interactions(parent_length)% &
           left_length_start(first_length-1)+1, &
           this%interactions(parent_length)%left_length_start(first_length)
        number = number+kinds
        number_of_attachments = number_of_attachments+kinds*( &
             this%interactions(parent_length)%attachment_start(pair)- &
             this%interactions(parent_length)%attachment_start(pair-1))
      end do
    end do
    if (number /= this%interactions(parent_length)%number_of_interactions .or. &
         number_of_attachments /= &
         this%interactions(parent_length)%number_of_attachments) &
         call fail('trace current DAG expanded schedule count mismatch')
  end subroutine build_interaction_level

  subroutine count_interaction_level(this, parent_length, word_level, &
       number, number_of_attachments, number_of_pairs, &
       number_of_pair_attachments)
    class(trace_current_dag_t), intent(in) :: this
    integer, intent(in) :: parent_length
    type(trace_word_level_t), intent(in) :: word_level(:)
    integer, intent(out) :: number, number_of_attachments, number_of_pairs
    integer, intent(out) :: number_of_pair_attachments

    integer(i64) :: attachments_64, interactions_64, pair_attachments_64
    integer(i64) :: pairs, pairs_64
    integer :: first, first_length, kinds, second, second_length, targets

    interactions_64 = 0_i64
    attachments_64 = 0_i64
    pairs_64 = 0_i64
    pair_attachments_64 = 0_i64
    if (.not. this%fixed_first_value) then
      do first_length = 1, parent_length-1
        second_length = parent_length-first_length
        pairs = full_disjoint_pair_count(this%degree_value, first_length, &
             second_length)
        kinds = number_of_interaction_kinds(first_length, second_length, &
             parent_length, this%degree_value)
        targets = full_interaction_target_count(first_length, second_length)
        interactions_64 = interactions_64+pairs*int(kinds, i64)
        attachments_64 = attachments_64+ &
             pairs*int(kinds*targets, i64)
        pairs_64 = pairs_64+pairs
        pair_attachments_64 = pair_attachments_64+pairs*int(targets, i64)
      end do
    else
      do first_length = 1, parent_length-1
        second_length = parent_length-first_length
        kinds = number_of_interaction_kinds(first_length, second_length, &
             parent_length, this%degree_value)
        do first = 1, size(word_level(first_length)%word, 2)
          do second = 1, size(word_level(second_length)%word, 2)
            if (iand(this%level(first_length)%subset_mask(first), &
                 this%level(second_length)%subset_mask(second)) /= 0) cycle
            if (word_level(first_length)%maximum(first) >= &
                 word_level(second_length)%maximum(second)) cycle
            targets = count_interaction_targets( &
                 word_level(first_length)%word(:, first), &
                 word_level(second_length)%word(:, second), .true.)
            interactions_64 = interactions_64+int(kinds, i64)
            attachments_64 = attachments_64+int(kinds*targets, i64)
            pairs_64 = pairs_64+1_i64
            pair_attachments_64 = pair_attachments_64+int(targets, i64)
          end do
        end do
      end do
    end if
    if (interactions_64 > int(huge(0), i64) .or. &
         attachments_64 > int(huge(0), i64) .or. &
         pairs_64 > int(huge(0), i64) .or. &
         pair_attachments_64 > int(huge(0), i64)) &
         call fail('trace current DAG interaction level exceeds array capacity')
    number = int(interactions_64)
    number_of_attachments = int(attachments_64)
    number_of_pairs = int(pairs_64)
    number_of_pair_attachments = int(pair_attachments_64)
  end subroutine count_interaction_level

  integer(i64) function full_disjoint_pair_count(degree, first_length, &
       second_length) result(number)
    integer, intent(in) :: degree, first_length, second_length

    integer :: parent_length, reversal_factor

    parent_length = first_length+second_length
    reversal_factor = 1
    if (first_length >= 2) reversal_factor = 2*reversal_factor
    if (second_length >= 2) reversal_factor = 2*reversal_factor
    ! In the full basis, P(degree,parent_length) orders the union labels.
    ! The largest label lies in the second child for the fraction
    ! second_length/parent_length; each nonsingleton child then has two
    ! orientations represented by one canonical word.
    number = falling_factorial_i64(degree, parent_length)* &
         int(second_length, i64)/int(parent_length*reversal_factor, i64)
  end function full_disjoint_pair_count

  integer function full_interaction_target_count(first_length, second_length) &
       result(number)
    integer, intent(in) :: first_length, second_length

    number = 1
    if (first_length >= 2) number = 2*number
    if (second_length >= 2) number = 2*number
  end function full_interaction_target_count

  integer function number_of_interaction_kinds(first_length, second_length, &
       parent_length, degree) result(number)
    integer, intent(in) :: first_length, second_length, parent_length, degree
    integer :: kind

    number = 0
    do kind = interaction_three_gluon, interaction_gluon_tensor
      if (interaction_allowed(kind, first_length, second_length, &
           parent_length, degree)) number = number+1
    end do
  end function number_of_interaction_kinds

  subroutine set_interaction_pair(level, pair, first, second)
    type(trace_interaction_level_t), intent(inout) :: level
    integer, intent(in) :: pair, first, second

    level%left_current(pair) = first
    level%right_current(pair) = second
  end subroutine set_interaction_pair

  subroutine append_interaction_targets(level, first, &
       second, parent_map, degree, rank_weight, cursor)
    type(trace_interaction_level_t), intent(inout) :: level
    integer, intent(in) :: first(:), second(:)
    integer, intent(in) :: parent_map(:), degree
    integer(i64), intent(in) :: rank_weight(:)
    integer, intent(inout) :: cursor

    integer :: first_switch, orientation, second_switch
    integer :: invert_first, invert_second, parent_rank, signed_target
    integer :: first_word(size(first)), parent(size(first)+size(second))
    integer :: second_word(size(second))
    do first_switch = 0, merge(1, 0, size(first) >= 2)
      do second_switch = 0, merge(1, 0, size(second) >= 2)
        do orientation = 0, 1
          invert_first = first_switch
          invert_second = second_switch
          if (invert_first == 0) then
            first_word = first
          else
            first_word = first(size(first):1:-1)
          end if
          if (invert_second == 0) then
            second_word = second
          else
            second_word = second(size(second):1:-1)
          end if
          if (orientation == 0) then
            parent = [first_word, second_word]
          else
            parent = [second_word, first_word]
          end if
          if (.not. canonical_word(parent)) cycle
          parent_rank = partial_permutation_rank(parent, degree, rank_weight)
          signed_target = parent_map(parent_rank)
          if (signed_target == 0) cycle
          if (signed_target <= 0) &
               call fail('canonical trace-current target has a negative map')
          if (interaction_target_is_negative(orientation, &
               first_switch, second_switch, size(first), size(second))) &
               signed_target = -signed_target
          cursor = cursor+1
          level%attachment_target(cursor) = signed_target
        end do
      end do
    end do
  end subroutine append_interaction_targets

  logical function interaction_target_is_negative(orientation, &
       first_switch, second_switch, first_length, second_length) result(negative)
    integer, intent(in) :: orientation, first_switch, second_switch
    integer, intent(in) :: first_length, second_length

    logical :: vertex_sign

    vertex_sign = orientation == 1
    if (first_switch == 1 .and. modulo(first_length, 2) == 0) &
         vertex_sign = .not. vertex_sign
    if (second_switch == 1 .and. modulo(second_length, 2) == 0) &
         vertex_sign = .not. vertex_sign

    ! Swapping the children changes the sign for every allowed interaction.
    negative = vertex_sign
  end function interaction_target_is_negative

  integer function count_interaction_targets(first, second, fixed_first) &
       result(number)
    integer, intent(in) :: first(:), second(:)
    logical, intent(in) :: fixed_first

    integer :: first_switch, orientation, second_switch
    integer :: invert_first, invert_second
    integer :: first_word(size(first)), parent(size(first)+size(second))
    integer :: second_word(size(second))

    if (.not. fixed_first) then
      number = 2
      if (size(first) >= 2) number = 2*number
      if (size(second) >= 2) number = 2*number
      number = number/2
      return
    end if

    number = 0
    do first_switch = 0, merge(1, 0, size(first) >= 2)
      do second_switch = 0, merge(1, 0, size(second) >= 2)
        do orientation = 0, 1
          invert_first = first_switch
          invert_second = second_switch
          if (invert_first == 0) then
            first_word = first
          else
            first_word = first(size(first):1:-1)
          end if
          if (invert_second == 0) then
            second_word = second
          else
            second_word = second(size(second):1:-1)
          end if
          if (orientation == 0) then
            parent = [first_word, second_word]
          else
            parent = [second_word, first_word]
          end if
          if (.not. retained_canonical_word(parent, fixed_first)) cycle
          number = number+1
        end do
      end do
    end do
  end function count_interaction_targets

  logical function interaction_allowed(kind, first_length, second_length, &
       parent_length, degree) result(allowed)
    integer, intent(in) :: kind, first_length, second_length
    integer, intent(in) :: parent_length, degree

    select case(kind)
    case(interaction_three_gluon)
      allowed = .true.
    case(interaction_tensor)
      allowed = parent_length < degree
    case(interaction_tensor_gluon)
      allowed = first_length >= 2
    case(interaction_gluon_tensor)
      allowed = second_length >= 2
    case default
      allowed = .false.
    end select
  end function interaction_allowed

  logical function canonical_word(word) result(canonical)
    integer, intent(in) :: word(:)

    integer :: maximum_position, minimum_position, position

    if (size(word) < 2) then
      canonical = .true.
    else
      minimum_position = 1
      maximum_position = 1
      do position = 2, size(word)
        if (word(position) < word(minimum_position)) &
             minimum_position = position
        if (word(position) > word(maximum_position)) &
             maximum_position = position
      end do
      canonical = minimum_position < maximum_position
    end if
  end function canonical_word

  logical function retained_canonical_word(word, fixed_first) result(retained)
    integer, intent(in) :: word(:)
    logical, intent(in) :: fixed_first

    retained = canonical_word(word)
    if (.not. retained .or. .not. fixed_first) return
    retained = .not. any(word == 1) .or. word(1) == 1
  end function retained_canonical_word

  subroutine evaluate_trace_current_dag(this, outgoing_momenta, &
       external_wavefunctions, terminal_wavefunction, coupling, amplitudes)
    class(trace_current_dag_t), intent(inout) :: this
    real(dp), intent(in) :: outgoing_momenta(0:, :)
    complex(dp), intent(in) :: external_wavefunctions(:, :)
    complex(dp), intent(in) :: terminal_wavefunction(:)
    real(dp), intent(in) :: coupling
    complex(dp), intent(out) :: amplitudes(:)

    complex(dp) :: tensor_value(6), terminal_value
    complex(dp) :: vector_term(4), vector_value(4)
    integer :: attachment, current, first, first_length, first_mask, length
    integer :: mask, order, pair, second, second_length, second_mask
    integer :: signed_target, target
    real(dp) :: amplitude_scale

    if (.not. this%initialized) call fail('trace current DAG is not initialized')
    if (size(outgoing_momenta, 1) /= 4 .or. &
         size(outgoing_momenta, 2) < this%degree_value) &
         call fail('wrong momentum array for trace current DAG')
    if (size(external_wavefunctions, 1) /= 4 .or. &
         size(external_wavefunctions, 2) < this%degree_value) &
         call fail('wrong wavefunction array for trace current DAG')
    if (size(terminal_wavefunction) /= 4) &
         call fail('wrong terminal wavefunction for trace current DAG')
    if (size(amplitudes) /= this%number_of_orders_value) &
         call fail('wrong amplitude array for trace current DAG')

    amplitude_scale = coupling**(this%degree_value-1)
    if (.not. ieee_is_finite(amplitude_scale)) &
         call fail('trace current DAG coupling power is not finite')
    call fill_subset_momenta(this, outgoing_momenta)
    do current = 1, this%degree_value
      this%level(1)%vector(:, current) = external_wavefunctions(:, current)
    end do
    do length = 2, this%degree_value
      this%level(length)%vector = cmplx(0.0_dp, 0.0_dp, dp)
      if (allocated(this%level(length)%tensor)) &
           this%level(length)%tensor = cmplx(0.0_dp, 0.0_dp, dp)
      if (length == this%degree_value) then
        do first_length = 1, length-1
          second_length = length-first_length
          do pair = this%interactions(length)% &
               left_length_start(first_length-1)+1, &
               this%interactions(length)%left_length_start(first_length)
          first = this%interactions(length)%left_current(pair)
          second = this%interactions(length)%right_current(pair)
          first_mask = this%level(first_length)%subset_mask(first)
          second_mask = this%level(second_length)%subset_mask(second)

          terminal_value = three_gluon_terminal( &
               this%level(first_length)%vector(:, first), &
               this%subset_momentum(:, first_mask), &
               this%level(second_length)%vector(:, second), &
               this%subset_momentum(:, second_mask), terminal_wavefunction)
          if (first_length >= 2) then
            terminal_value = terminal_value+aux_tensor_gluon_terminal( &
                 this%level(first_length)%tensor(:, first), &
                 this%level(second_length)%vector(:, second), &
                 terminal_wavefunction)
          end if
          if (second_length >= 2) then
            terminal_value = terminal_value+gluon_aux_tensor_terminal( &
                 this%level(first_length)%vector(:, first), &
                 this%level(second_length)%tensor(:, second), &
                 terminal_wavefunction)
          end if
          do attachment = this%interactions(length)% &
               attachment_start(pair-1)+1, &
               this%interactions(length)%attachment_start(pair)
            signed_target = this%interactions(length)% &
                 attachment_target(attachment)
            target = abs(signed_target)
            if (signed_target > 0) then
              this%level(length)%vector(1, target) = &
                   this%level(length)%vector(1, target)+terminal_value
            else
              this%level(length)%vector(1, target) = &
                   this%level(length)%vector(1, target)-terminal_value
            end if
          end do
          end do
        end do
      else
        do first_length = 1, length-1
          second_length = length-first_length
          do pair = this%interactions(length)% &
               left_length_start(first_length-1)+1, &
               this%interactions(length)%left_length_start(first_length)
          first = this%interactions(length)%left_current(pair)
          second = this%interactions(length)%right_current(pair)
          first_mask = this%level(first_length)%subset_mask(first)
          second_mask = this%level(second_length)%subset_mask(second)

          call three_gluon(this%level(first_length)%vector(:, first), &
               this%subset_momentum(:, first_mask), &
               this%level(second_length)%vector(:, second), &
               this%subset_momentum(:, second_mask), vector_value)
          if (first_length >= 2) then
            call aux_tensor_gluon_to_gluon( &
                 this%level(first_length)%tensor(:, first), &
                 this%level(second_length)%vector(:, second), vector_term)
            vector_value = vector_value+vector_term
          end if
          if (second_length >= 2) then
            call gluon_aux_tensor_to_gluon( &
                 this%level(first_length)%vector(:, first), &
                 this%level(second_length)%tensor(:, second), vector_term)
            vector_value = vector_value+vector_term
          end if
          call two_gluon_to_aux_tensor( &
               this%level(first_length)%vector(:, first), &
               this%level(second_length)%vector(:, second), tensor_value)
          do attachment = this%interactions(length)% &
               attachment_start(pair-1)+1, &
               this%interactions(length)%attachment_start(pair)
            signed_target = this%interactions(length)% &
                 attachment_target(attachment)
            target = abs(signed_target)
            if (signed_target > 0) then
              this%level(length)%vector(:, target) = &
                   this%level(length)%vector(:, target)+vector_value
              this%level(length)%tensor(:, target) = &
                   this%level(length)%tensor(:, target)+tensor_value
            else
              this%level(length)%vector(:, target) = &
                   this%level(length)%vector(:, target)-vector_value
              this%level(length)%tensor(:, target) = &
                   this%level(length)%tensor(:, target)-tensor_value
            end if
          end do
          end do
        end do

        do current = 1, this%level(length)%number_of_currents
          mask = this%level(length)%subset_mask(current)
          call massless_vector_propagator( &
               this%level(length)%vector(:, current), &
               this%subset_momentum(:, mask))
        end do
      end if
    end do

    do order = 1, size(amplitudes)
      signed_target = this%signed_order_current(order)
      terminal_value = this%level(this%degree_value)% &
           vector(1, abs(signed_target))
      if (signed_target < 0) terminal_value = -terminal_value
      amplitudes(order) = terminal_value
    end do
    if (abs(amplitude_scale-1.0_dp) > 0.0_dp) &
         amplitudes = amplitude_scale*amplitudes
  end subroutine evaluate_trace_current_dag

  subroutine fill_subset_momenta(this, outgoing_momenta)
    class(trace_current_dag_t), intent(inout) :: this
    real(dp), intent(in) :: outgoing_momenta(0:, :)

    integer :: leg, mask, previous_mask

    this%subset_momentum(:, 0) = 0.0_dp
    do mask = 1, ubound(this%subset_momentum, 2)
      leg = trailz(mask)+1
      previous_mask = ibclr(mask, leg-1)
      this%subset_momentum(:, mask) = &
           this%subset_momentum(:, previous_mask)+outgoing_momenta(:, leg)
    end do
  end subroutine fill_subset_momenta

  subroutine build_order_map(this, words)
    class(trace_current_dag_t), intent(inout) :: this
    integer, intent(in) :: words(:, :)

    integer, allocatable :: signed_map(:), tail(:), word(:)
    integer(i64), allocatable :: rank_weight(:)
    integer :: order

    call build_signed_word_map(this%degree_value, this%degree_value, &
         this%fixed_first_value, words, signed_map)
    allocate(this%signed_order_current(this%number_of_orders_value))
    allocate(word(this%degree_value))
    allocate(rank_weight(this%degree_value))
    call fill_partial_permutation_weights(this%degree_value, &
         this%degree_value, rank_weight)
    if (this%fixed_first_value) then
      allocate(tail(this%degree_value-1))
      word(1) = 1
      tail = [(order, order=2, this%degree_value)]
      do order = 1, this%number_of_orders_value
        word(2:) = tail
        this%signed_order_current(order) = signed_map( &
             partial_permutation_rank(word, this%degree_value, rank_weight))
        if (this%signed_order_current(order) <= 0) &
             call fail('fixed-first trace current DAG order was not mapped')
        if (order < this%number_of_orders_value) call next_permutation(tail)
      end do
    else
      word = [(order, order=1, this%degree_value)]
      do order = 1, this%number_of_orders_value
        this%signed_order_current(order) = signed_map( &
             partial_permutation_rank(word, this%degree_value, rank_weight))
        if (this%signed_order_current(order) == 0) &
             call fail('trace current DAG order was not mapped')
        if (order < this%number_of_orders_value) call next_permutation(word)
      end do
    end if
  end subroutine build_order_map

  subroutine build_signed_word_map(degree, length, fixed_first, words, &
       signed_map)
    integer, intent(in) :: degree, length
    logical, intent(in) :: fixed_first
    integer, intent(in) :: words(:, :)
    integer, allocatable, intent(out) :: signed_map(:)

    integer, allocatable :: reverse_word(:)
    integer :: current, expected_mapped, rank, reverse_rank, reverse_sign
    integer(i64), allocatable :: rank_weight(:)

    allocate(signed_map(int(falling_factorial_i64(degree, length))))
    signed_map = 0
    allocate(reverse_word(length))
    allocate(rank_weight(length))
    call fill_partial_permutation_weights(degree, length, rank_weight)
    reverse_sign = merge(1, -1, modulo(length-1, 2) == 0)
    do current = 1, size(words, 2)
      rank = partial_permutation_rank(words(:, current), degree, rank_weight)
      reverse_word = words(length:1:-1, current)
      reverse_rank = partial_permutation_rank(reverse_word, degree, rank_weight)
      signed_map(rank) = current
      signed_map(reverse_rank) = reverse_sign*current
    end do
    if (fixed_first) then
      expected_mapped = size(words, 2)
      if (length >= 2) expected_mapped = 2*expected_mapped
      if (count(signed_map /= 0) /= expected_mapped) &
           call fail('incomplete fixed-first trace-current signed-word map')
    else if (any(signed_map == 0)) then
      call fail('incomplete trace current DAG signed-word map')
    end if
  end subroutine build_signed_word_map

  subroutine enumerate_canonical_words(degree, length, fixed_first, words)
    integer, intent(in) :: degree, length
    logical, intent(in) :: fixed_first
    integer, allocatable, intent(out) :: words(:, :)

    integer, allocatable :: word(:)
    logical, allocatable :: used(:)
    integer :: count, cursor

    if (length < 1 .or. length > degree) &
         call fail('invalid trace-current word length')
    if (fixed_first .and. length >= 2) then
      count = int(falling_factorial_i64(degree-1, length-1))
      if (length <= degree-1) count = count+ &
           int(falling_factorial_i64(degree-1, length)/2_i64)
    else
      count = int(falling_factorial_i64(degree, length))
      if (length >= 2) count = count/2
    end if
    allocate(words(length, count), word(length), used(degree))
    used = .false.
    cursor = 0
    call append_words(1)
    if (cursor /= count) call fail('trace current DAG word count mismatch')
  contains
    recursive subroutine append_words(position)
      integer, intent(in) :: position
      integer :: label

      if (position > length) then
        if (.not. canonical_word(word)) return
        cursor = cursor+1
        words(:, cursor) = word
        return
      end if
      do label = 1, degree
        if (used(label)) cycle
        if (fixed_first .and. position > 1 .and. label == 1) then
          if (word(1) /= 1) cycle
        end if
        word(position) = label
        used(label) = .true.
        call append_words(position+1)
        used(label) = .false.
      end do
    end subroutine append_words
  end subroutine enumerate_canonical_words

  subroutine fill_subset_masks(words, masks)
    integer, intent(in) :: words(:, :)
    integer, intent(out) :: masks(:)

    integer :: current, position

    if (size(masks) /= size(words, 2)) &
         call fail('trace-current subset-mask dimension mismatch')
    masks = 0
    do current = 1, size(words, 2)
      do position = 1, size(words, 1)
        masks(current) = ibset(masks(current), words(position, current)-1)
      end do
    end do
  end subroutine fill_subset_masks

  subroutine fill_partial_permutation_weights(degree, length, weight)
    integer, intent(in) :: degree, length
    integer(i64), intent(out) :: weight(:)

    integer :: position

    if (size(weight) /= length) &
         call fail('trace-current rank-weight dimension mismatch')
    do position = 1, length
      weight(position) = falling_factorial_i64(degree-position, &
           length-position)
    end do
  end subroutine fill_partial_permutation_weights

  integer function partial_permutation_rank(word, degree, weight) result(rank)
    integer, intent(in) :: word(:), degree
    integer(i64), intent(in) :: weight(:)

    integer :: label, lower_labels, position, smaller, used_mask
    integer(i64) :: rank_64

    used_mask = 0
    rank_64 = 1_i64
    if (size(weight) /= size(word)) &
         call fail('trace-current rank-weight dimension mismatch')
    do position = 1, size(word)
      label = word(position)
      if (label < 1 .or. label > degree) &
           call fail('invalid trace-current word')
      if (btest(used_mask, label-1)) call fail('invalid trace-current word')
      lower_labels = shiftl(1, label-1)-1
      smaller = label-1-popcnt(iand(used_mask, lower_labels))
      rank_64 = rank_64+int(smaller, i64)*weight(position)
      used_mask = ibset(used_mask, label-1)
    end do
    if (rank_64 > int(huge(0), i64)) &
         call fail('trace-current word rank exceeds array capacity')
    rank = int(rank_64)
  end function partial_permutation_rank

  integer(i64) function falling_factorial_i64(number, length) result(value)
    integer, intent(in) :: number, length

    integer :: factor

    if (number < 0 .or. length < 0 .or. length > number) &
         call fail('invalid falling factorial in trace current DAG')
    value = 1_i64
    do factor = number-length+1, number
      if (value > huge(value)/int(factor, i64)) &
           call fail('trace current DAG factorial overflows 64 bits')
      value = value*int(factor, i64)
    end do
  end function falling_factorial_i64

  subroutine next_permutation(permutation)
    integer, intent(inout) :: permutation(:)

    integer :: left, right, temporary

    left = size(permutation)-1
    do while (left >= 1)
      if (permutation(left) < permutation(left+1)) exit
      left = left-1
    end do
    if (left < 1) call fail('trace-current permutation wrapped unexpectedly')
    right = size(permutation)
    do while (permutation(right) < permutation(left))
      right = right-1
    end do
    temporary = permutation(left)
    permutation(left) = permutation(right)
    permutation(right) = temporary
    permutation(left+1:) = permutation(size(permutation):left+1:-1)
  end subroutine next_permutation

  integer function number_of_word_slots(this) result(number)
    class(trace_current_dag_t), intent(in) :: this
    integer :: length

    if (.not. this%initialized) call fail('trace current DAG is not initialized')
    number = 0
    do length = 1, this%degree_value
      number = number+this%level(length)%number_of_currents
    end do
  end function number_of_word_slots

  integer function number_of_interactions(this) result(number)
    class(trace_current_dag_t), intent(in) :: this
    integer :: length

    if (.not. this%initialized) call fail('trace current DAG is not initialized')
    number = 0
    do length = 2, this%degree_value
      number = number+this%interactions(length)%number_of_pairs
    end do
  end function number_of_interactions

  integer function number_of_attachments(this) result(number)
    class(trace_current_dag_t), intent(in) :: this
    integer :: length

    if (.not. this%initialized) call fail('trace current DAG is not initialized')
    number = 0
    do length = 2, this%degree_value
      number = number+this%interactions(length)%number_of_pair_attachments
    end do
  end function number_of_attachments

  integer function number_of_orders(this) result(number)
    class(trace_current_dag_t), intent(in) :: this

    if (.not. this%initialized) call fail('trace current DAG is not initialized')
    number = this%number_of_orders_value
  end function number_of_orders

  subroutine ensure_schedule_fits_default_integer(degree, fixed_first)
    integer, intent(in) :: degree
    logical, intent(in) :: fixed_first

    integer(i64) :: attachments, currents
    integer :: length, terms_per_current

    attachments = 0_i64
    do length = 2, degree
      if (fixed_first) then
        currents = falling_factorial_i64(degree-1, length-1)
        if (length <= degree-1) currents = currents+ &
             falling_factorial_i64(degree-1, length)/2_i64
      else
        currents = falling_factorial_i64(degree, length)/2_i64
      end if
      terms_per_current = 3*length-5
      if (length < degree) terms_per_current = terms_per_current+length-1
      if (currents > (int(huge(0), i64)-attachments)/ &
           int(terms_per_current, i64)) &
           call fail('trace current DAG attachment count exceeds array capacity')
      attachments = attachments+currents*int(terms_per_current, i64)
    end do
  end subroutine ensure_schedule_fits_default_integer

end module trace_current_dag
