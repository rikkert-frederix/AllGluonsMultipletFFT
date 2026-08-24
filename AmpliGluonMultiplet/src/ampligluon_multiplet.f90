module ampligluon_multiplet
  use iso_c_binding, only: c_bool
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ampligluon_multiplet_kinds, only: dp, fail, i8
  use multiplet_paths, only: path_catalog_t
  use multiplet_mhv, only: multiplet_mhv_t
  use qcd_kinematics, only: aux_tensor_gluon_to_gluon, &
       external_massless_vector, gluon_aux_tensor_to_gluon, &
       lorentz_contract, massless_vector_propagator, &
       massless_vector_propagator_copy, three_gluon, two_gluon_to_aux_tensor
  use recoupling_plan, only: apply_sparse_add, apply_sparse_pair_add, &
       apply_sparse_scalar_add, recoupling_system_t
  use wigner_table, only: wigner_table_t
  implicit none
  private

  real(dp), parameter :: adjoint_dimension = 8.0_dp
  real(dp), parameter :: normalized_bracket_factor = sqrt(6.0_dp)
  real(dp), parameter :: momentum_tolerance = 2.0e-9_dp

  type :: subset_current_t
    integer :: length = 0
    real(dp) :: momentum(0:3) = 0.0_dp
    complex(dp), allocatable :: gluon(:, :)
    complex(dp), allocatable :: tensor(:, :)
  end type subset_current_t

  type :: subset_evaluation_t
    integer(i8), allocatable :: left_pattern(:)
    integer(i8), allocatable :: right_pattern(:)
  end type subset_evaluation_t

  type :: recoupling_workspace_t
    complex(dp), allocatable :: gluon(:, :, :, :)
    complex(dp), allocatable :: tensor(:, :, :, :)
    complex(dp), allocatable :: scalar(:, :, :)
    logical(c_bool), allocatable :: initialized(:, :, :)
  end type recoupling_workspace_t

  type, public :: ampligluon_multiplet_t
    private
    logical :: initialized = .false.
    logical :: generic_initialized = .false.
    integer :: total_gluons = 0
    integer :: current_gluons = 0
    integer(i8) :: maximum_pattern = 0_i8
    type(wigner_table_t) :: table
    type(path_catalog_t) :: catalog
    type(recoupling_system_t) :: recoupling
    type(multiplet_mhv_t) :: mhv
    type(subset_current_t), allocatable :: currents(:)
    type(subset_evaluation_t), allocatable :: evaluations(:)
    type(recoupling_workspace_t), allocatable :: workspaces(:)
  contains
    procedure, public :: evaluate => evaluate_ampligluon_multiplet
    procedure, public :: initialize => initialize_ampligluon_multiplet
    procedure, public :: number_of_basis_amplitudes
    procedure, public :: number_of_final_gluons
  end type ampligluon_multiplet_t

contains

  subroutine initialize_ampligluon_multiplet(this, final_gluons, table_file)
    class(ampligluon_multiplet_t), intent(inout) :: this
    integer, intent(in) :: final_gluons
    character(len=*), intent(in) :: table_file

    integer :: required_prefix_gluons

    if (final_gluons < 2) call fail('gg scattering needs at least two final gluons')
    this%total_gluons = final_gluons+2
    this%current_gluons = this%total_gluons-1
    if (this%current_gluons >= bit_size(this%maximum_pattern)-1) &
         call fail('requested process exceeds the bit-mask capacity')
    if (this%current_gluons >= bit_size(0)-1) &
         call fail('requested process exceeds the array-index capacity')

    call this%table%load(table_file)
    ! Closing the current with its final adjoint makes a chain of
    ! total_gluons adjoints.  Two-sided reachability of every cut gives the
    ! sufficient prefix cutoff used by this one-sided path construction.
    required_prefix_gluons = this%total_gluons/2
    if (required_prefix_gluons > this%table%max_prefix_gluons) then
      call fail('Wigner table is too shallow for the requested process')
    end if
    call this%catalog%build(this%table, this%current_gluons)
    call this%recoupling%build_swaps(this%table, this%catalog)
    call this%mhv%initialize(this%catalog)

    this%maximum_pattern = shiftl(1_i8, this%current_gluons)-1_i8
    if (allocated(this%currents)) deallocate(this%currents)
    if (allocated(this%evaluations)) deallocate(this%evaluations)
    if (allocated(this%workspaces)) deallocate(this%workspaces)
    this%generic_initialized = .false.
    this%initialized = .true.
  end subroutine initialize_ampligluon_multiplet

  subroutine initialize_generic_workspaces(this)
    class(ampligluon_multiplet_t), intent(inout) :: this

    integer :: length, number_of_paths, pattern, workspace_slots

    if (this%generic_initialized) return
    call this%recoupling%build_plans(this%catalog)
    allocate(this%currents(int(this%maximum_pattern)))
    allocate(this%evaluations(int(this%maximum_pattern)))
    allocate(this%workspaces(this%current_gluons))

    do pattern = 1, int(this%maximum_pattern)
      length = popcnt(pattern)
      number_of_paths = this%catalog%spaces(length)%number_of_paths
      this%currents(pattern)%length = length
      if (length < this%current_gluons) then
        allocate(this%currents(pattern)%gluon(4, number_of_paths))
        allocate(this%currents(pattern)%tensor(6, number_of_paths))
        this%currents(pattern)%gluon = (0.0_dp, 0.0_dp)
        this%currents(pattern)%tensor = (0.0_dp, 0.0_dp)
      end if
      if (length >= 2) call build_subset_evaluation(this, pattern)
    end do

    do length = 2, this%current_gluons
      number_of_paths = this%catalog%spaces(length)%number_of_paths
      workspace_slots = this%recoupling%plans(length)%max_layer_width
      allocate(this%workspaces(length)%initialized(number_of_paths, &
                                                   workspace_slots, 2))
      if (length < this%current_gluons) then
        allocate(this%workspaces(length)%gluon(4, number_of_paths, &
                                               workspace_slots, 2))
        allocate(this%workspaces(length)%tensor(6, number_of_paths, &
                                                workspace_slots, 2))
      else
        allocate(this%workspaces(length)%scalar(number_of_paths, &
                                                workspace_slots, 2))
      end if
    end do
    this%generic_initialized = .true.
  end subroutine initialize_generic_workspaces

  subroutine build_subset_evaluation(this, pattern)
    class(ampligluon_multiplet_t), intent(inout) :: this
    integer, intent(in) :: pattern

    integer(i8) :: actual_left, relative_left
    integer :: external_leg, index_partition, rank

    associate(evaluation => this%evaluations(pattern), &
              plan => this%recoupling%plans(popcnt(pattern)))
      allocate(evaluation%left_pattern(plan%number_of_partitions))
      allocate(evaluation%right_pattern(plan%number_of_partitions))
      do index_partition = 1, plan%number_of_partitions
        relative_left = plan%partitions(index_partition)%left_pattern
        actual_left = 0_i8
        rank = 0
        do external_leg = 1, this%current_gluons
          if (.not. btest(pattern, external_leg-1)) cycle
          rank = rank+1
          if (btest(relative_left, rank-1)) &
               actual_left = ibset(actual_left, external_leg-1)
        end do
        evaluation%left_pattern(index_partition) = actual_left
        evaluation%right_pattern(index_partition) = &
             ieor(int(pattern, i8), actual_left)
      end do
    end associate
  end subroutine build_subset_evaluation

  subroutine evaluate_ampligluon_multiplet(this, momenta, helicities, &
                                           basis_amplitudes, matrix_element_squared, &
                                           strong_coupling, average_initial_colours, &
                                           use_mhv_optimization)
    class(ampligluon_multiplet_t), intent(inout) :: this
    real(dp), intent(in) :: momenta(0:, :)
    integer, intent(in) :: helicities(:)
    complex(dp), allocatable, intent(inout) :: basis_amplitudes(:)
    real(dp), intent(out) :: matrix_element_squared
    real(dp), intent(in), optional :: strong_coupling
    logical, intent(in), optional :: average_initial_colours
    logical, intent(in), optional :: use_mhv_optimization

    complex(dp) :: closing_wavefunction(4)
    real(dp) :: outgoing_momenta(0:3, size(momenta, 2))
    real(dp) :: coupling
    integer(i8) :: external_pattern
    integer :: external_leg, length, outgoing_positive, pattern, remainder
    integer :: top_pattern
    logical :: average_colours, enable_mhv_optimization, used_mhv

    if (.not. this%initialized) call fail('amplitude object is not initialized')
    call validate_event(this, momenta, helicities)
    coupling = 1.0_dp
    if (present(strong_coupling)) coupling = strong_coupling
    if (.not. ieee_is_finite(coupling)) call fail('strong coupling is not finite')
    average_colours = .false.
    if (present(average_initial_colours)) average_colours = average_initial_colours
    enable_mhv_optimization = .true.
    if (present(use_mhv_optimization)) &
         enable_mhv_optimization = use_mhv_optimization

    call ensure_basis_output(this, basis_amplitudes)

    ! Every connected pure-gluon tree contains exactly N-2 powers of the
    ! strong coupling. Handle zero (and numerically irrelevant subnormal)
    ! coupling before either recursion; this also avoids forming harmless
    ! 0/0 analytic ratios at singular kinematics.
    if (abs(coupling) < tiny(1.0_dp)) then
      basis_amplitudes = (0.0_dp, 0.0_dp)
      matrix_element_squared = 0.0_dp
      return
    end if

    ! Pure-gluon tree amplitudes vanish when their all-outgoing helicities
    ! contain fewer than two states of either sign.  Crossing flips the two
    ! incoming helicities.  Return the exact zero before constructing any
    ! Berends--Giele currents, while still providing the documented basis
    ! vector and performing all event/coupling validation above.
    outgoing_positive = count(helicities(1:2) < 0) + &
                        count(helicities(3:this%total_gluons) > 0)
    if (outgoing_positive < 2 .or. &
        outgoing_positive > this%total_gluons-2) then
      basis_amplitudes = (0.0_dp, 0.0_dp)
      matrix_element_squared = 0.0_dp
      return
    end if

    outgoing_momenta(:, 1:2) = -momenta(:, 1:2)
    outgoing_momenta(:, 3:this%total_gluons) = &
         momenta(:, 3:this%total_gluons)

    ! The MHV and anti-MHV sectors admit an exact inverse-soft recursion
    ! directly in this same orthonormal multiplet basis.  Its colour maps are
    ! Wigner-swap/antisymmetric-vertex sandwiches; if a spinor denominator is
    ! unsafe, the general multiplet Berends--Giele recursion remains the
    ! transparent fallback.
    used_mhv = .false.
    if (enable_mhv_optimization) &
         call this%mhv%evaluate(this%recoupling, outgoing_momenta, helicities, &
                                coupling, basis_amplitudes, used_mhv)
    if (used_mhv) then
      matrix_element_squared = &
           sum(real(basis_amplitudes*conjg(basis_amplitudes), dp))
      if (average_colours) matrix_element_squared = &
           matrix_element_squared/(adjoint_dimension**2)
      return
    end if

    call initialize_generic_workspaces(this)

    do pattern = 1, int(this%maximum_pattern)
      external_leg = trailz(pattern)+1
      remainder = ibclr(pattern, external_leg-1)
      this%currents(pattern)%momentum = outgoing_momenta(:, external_leg)
      if (remainder > 0) this%currents(pattern)%momentum = &
           this%currents(pattern)%momentum + this%currents(remainder)%momentum
    end do

    do external_leg = 1, this%current_gluons
      external_pattern = shiftl(1_i8, external_leg-1)
      call external_massless_vector(outgoing_momenta(:, external_leg), &
           helicities(external_leg), &
           this%currents(int(external_pattern))%gluon(:, 1))
      this%currents(int(external_pattern))%tensor = (0.0_dp, 0.0_dp)
    end do

    do length = 2, this%current_gluons-1
      do pattern = 1, int(this%maximum_pattern)
        if (this%currents(pattern)%length /= length) cycle
        call evaluate_subset(this, pattern)
      end do
    end do

    call external_massless_vector(outgoing_momenta(:, this%total_gluons), &
                                  helicities(this%total_gluons), &
                                  closing_wavefunction)
    top_pattern = int(this%maximum_pattern)
    call evaluate_terminal_subset(this, top_pattern, coupling, &
                                  closing_wavefunction, basis_amplitudes)
    matrix_element_squared = sum(real(basis_amplitudes*conjg(basis_amplitudes), dp))
    if (average_colours) matrix_element_squared = &
         matrix_element_squared/(adjoint_dimension**2)
  end subroutine evaluate_ampligluon_multiplet

  subroutine ensure_basis_output(this, basis_amplitudes)
    class(ampligluon_multiplet_t), intent(in) :: this
    complex(dp), allocatable, intent(inout) :: basis_amplitudes(:)

    integer :: required_size

    required_size = &
         this%catalog%spaces(this%current_gluons)%number_of_paths
    if (allocated(basis_amplitudes)) then
      if (size(basis_amplitudes) /= required_size) deallocate(basis_amplitudes)
    end if
    if (.not. allocated(basis_amplitudes)) allocate(basis_amplitudes(required_size))
  end subroutine ensure_basis_output

  subroutine evaluate_subset(this, pattern)
    class(ampligluon_multiplet_t), intent(inout) :: this
    integer, intent(in) :: pattern

    integer :: active_index, canonical_slot, current_buffer, inversion
    integer :: layer_entry, next_buffer, node, path, source_slot
    integer :: swap_position, target_node, target_slot, total_length
    logical :: build_tensor

    total_length = this%currents(pattern)%length
    build_tensor = total_length < this%current_gluons
    associate(workspace => this%workspaces(total_length), &
              plan => this%recoupling%plans(total_length))
      ! Sorting edges always join adjacent inversion layers.  Seed the
      ! highest layer, then alternate two compact layer buffers.  Seeds for
      ! the destination are accumulated before incoming edges, matching the
      ! addition order of the former all-node workspace.
      current_buffer = 1
      call clear_workspace_layer(plan%max_inversions, current_buffer)
      call accumulate_seed_layer(plan%max_inversions, current_buffer)

      do inversion = plan%max_inversions, 1, -1
        next_buffer = 3-current_buffer
        call clear_workspace_layer(inversion-1, next_buffer)
        call accumulate_seed_layer(inversion-1, next_buffer)

        do layer_entry = plan%layer_node_start(inversion), &
                         plan%layer_node_start(inversion+1)-1
          node = plan%layer_nodes(layer_entry)
          swap_position = plan%nodes(node)%swap_position
          target_node = plan%nodes(node)%next_node
          if (swap_position < 1 .or. target_node < 1) &
               call fail('incomplete recoupling ordering edge')
          source_slot = plan%nodes(node)%layer_slot
          target_slot = plan%nodes(target_node)%layer_slot
          if (build_tensor) then
            call apply_sparse_pair_add( &
                 this%recoupling%swaps(total_length)%matrices(swap_position), &
                 workspace%gluon(:, :, source_slot, current_buffer), &
                 workspace%tensor(:, :, source_slot, current_buffer), &
                 workspace%gluon(:, :, target_slot, next_buffer), &
                 workspace%tensor(:, :, target_slot, next_buffer), &
                 plan%nodes(node)%active_path, &
                 workspace%initialized(:, target_slot, next_buffer))
          else
            call apply_sparse_add( &
                 this%recoupling%swaps(total_length)%matrices(swap_position), &
                 workspace%gluon(:, :, source_slot, current_buffer), &
                 workspace%gluon(:, :, target_slot, next_buffer), &
                 plan%nodes(node)%active_path, &
                 workspace%initialized(:, target_slot, next_buffer))
          end if
        end do
        current_buffer = next_buffer
      end do

      canonical_slot = plan%nodes(plan%canonical_node)%layer_slot
      if (size(plan%nodes(plan%canonical_node)%active_path) == &
          plan%number_of_paths) then
        call massless_vector_propagator_copy( &
             workspace%gluon(:, :, canonical_slot, current_buffer), &
             this%currents(pattern)%gluon, this%currents(pattern)%momentum)
        this%currents(pattern)%tensor = &
             workspace%tensor(:, :, canonical_slot, current_buffer)
      else
        do active_index = 1, &
             size(plan%nodes(plan%canonical_node)%active_path)
          path = plan%nodes(plan%canonical_node)%active_path(active_index)
          this%currents(pattern)%gluon(:, path) = &
               workspace%gluon(:, path, canonical_slot, current_buffer)
          this%currents(pattern)%tensor(:, path) = &
               workspace%tensor(:, path, canonical_slot, current_buffer)
        end do
        call massless_vector_propagator(this%currents(pattern)%gluon, &
                                        this%currents(pattern)%momentum)
      end if
    end associate

  contains

    subroutine clear_workspace_layer(inversion_layer, buffer)
      integer, intent(in) :: inversion_layer, buffer

      integer :: active, entry_layer, path_index, path_slot, workspace_node

      associate(workspace => this%workspaces(total_length), &
                plan => this%recoupling%plans(total_length))
        do entry_layer = plan%layer_node_start(inversion_layer), &
                         plan%layer_node_start(inversion_layer+1)-1
          workspace_node = plan%layer_nodes(entry_layer)
          path_slot = plan%nodes(workspace_node)%layer_slot
          if (size(plan%nodes(workspace_node)%active_path) == &
              plan%number_of_paths) then
            workspace%initialized(:, path_slot, buffer) = .false.
          else
            do active = 1, size(plan%nodes(workspace_node)%active_path)
              path_index = plan%nodes(workspace_node)%active_path(active)
              workspace%initialized(path_index, path_slot, buffer) = .false.
            end do
          end if
        end do
      end associate
    end subroutine clear_workspace_layer

    subroutine accumulate_seed_layer(inversion_layer, buffer)
      integer, intent(in) :: inversion_layer, buffer

      complex(dp) :: vertex_gluon(4), vertex_piece(4), vertex_tensor(6)
      integer :: initial_node, initial_slot, index_partition, layer_seed
      integer :: left_active, left_index, left_path, left_size
      integer :: parent_path, right_active, right_index, right_path
      integer :: signed_parent_path

      associate(workspace => this%workspaces(total_length), &
                evaluation => this%evaluations(pattern), &
                plan => this%recoupling%plans(total_length))
        do layer_seed = plan%layer_partition_start(inversion_layer), &
                        plan%layer_partition_start(inversion_layer+1)-1
          index_partition = plan%layer_partitions(layer_seed)
          ! Only one orientation of each bipartition is stored (the lowest
          ! ranked leg is on the left).  The antisymmetry of both the colour
          ! bracket and the kinematic cubic/auxiliary vertices makes their
          ! product symmetric, so the opposite orientation would duplicate
          ! it.
          left_index = int(evaluation%left_pattern(index_partition))
          right_index = int(evaluation%right_pattern(index_partition))
          left_size = plan%partitions(index_partition)%left_size
          initial_node = plan%partitions(index_partition)%initial_node
          initial_slot = plan%nodes(initial_node)%layer_slot
          associate(seed => plan%seeds(left_size), &
                    left => this%currents(left_index), &
                    right => this%currents(right_index), &
                    left_support => &
                         this%recoupling%canonical_support(left_size), &
                    right_support => &
                         this%recoupling%canonical_support( &
                              total_length-left_size))
            do right_active = 1, size(right_support%active_path)
              right_path = right_support%active_path(right_active)
              do left_active = 1, size(left_support%active_path)
                left_path = left_support%active_path(left_active)
                signed_parent_path = &
                     seed%signed_parent_path(left_path, right_path)
                parent_path = abs(signed_parent_path)

                call three_gluon(left%gluon(:, left_path), left%momentum, &
                                 right%gluon(:, right_path), right%momentum, &
                                 vertex_gluon)

                if (left%length >= 2) then
                  call aux_tensor_gluon_to_gluon( &
                       left%tensor(:, left_path), &
                       right%gluon(:, right_path), vertex_piece)
                  vertex_gluon = vertex_gluon+vertex_piece
                end if
                if (right%length >= 2) then
                  call gluon_aux_tensor_to_gluon( &
                       left%gluon(:, left_path), &
                       right%tensor(:, right_path), vertex_piece)
                  vertex_gluon = vertex_gluon+vertex_piece
                end if
                if (.not. workspace%initialized( &
                    parent_path, initial_slot, buffer)) then
                  if (signed_parent_path > 0) then
                    workspace%gluon(:, parent_path, initial_slot, buffer) = &
                         vertex_gluon
                  else
                    workspace%gluon(:, parent_path, initial_slot, buffer) = &
                         -vertex_gluon
                  end if
                else
                  if (signed_parent_path > 0) then
                    workspace%gluon(:, parent_path, initial_slot, buffer) = &
                         workspace%gluon(:, parent_path, initial_slot, buffer) + &
                         vertex_gluon
                  else
                    workspace%gluon(:, parent_path, initial_slot, buffer) = &
                         workspace%gluon(:, parent_path, initial_slot, buffer) - &
                         vertex_gluon
                  end if
                end if

                if (build_tensor) then
                  call two_gluon_to_aux_tensor( &
                       left%gluon(:, left_path), right%gluon(:, right_path), &
                       vertex_tensor)
                  if (.not. workspace%initialized( &
                      parent_path, initial_slot, buffer)) then
                    if (signed_parent_path > 0) then
                      workspace%tensor(:, parent_path, initial_slot, buffer) = &
                           vertex_tensor
                    else
                      workspace%tensor(:, parent_path, initial_slot, buffer) = &
                           -vertex_tensor
                    end if
                  else
                    if (signed_parent_path > 0) then
                      workspace%tensor(:, parent_path, initial_slot, buffer) = &
                           workspace%tensor(:, parent_path, &
                                            initial_slot, buffer) + vertex_tensor
                    else
                      workspace%tensor(:, parent_path, initial_slot, buffer) = &
                           workspace%tensor(:, parent_path, &
                                            initial_slot, buffer) - vertex_tensor
                    end if
                  end if
                end if
                workspace%initialized(parent_path, initial_slot, buffer) = &
                     .true.
              end do
            end do
          end associate
        end do
      end associate
    end subroutine accumulate_seed_layer

  end subroutine evaluate_subset

  subroutine evaluate_terminal_subset(this, pattern, coupling, &
                                      closing_wavefunction, basis_amplitudes)
    class(ampligluon_multiplet_t), intent(inout) :: this
    integer, intent(in) :: pattern
    real(dp), intent(in) :: coupling
    complex(dp), intent(in) :: closing_wavefunction(4)
    complex(dp), intent(out) :: basis_amplitudes(:)

    real(dp) :: amplitude_scale
    integer :: active_index, canonical_slot, current_buffer, inversion
    integer :: layer_entry, next_buffer, node, path, source_slot
    integer :: swap_position, target_node, target_slot, total_length

    total_length = this%current_gluons
    amplitude_scale = sqrt(adjoint_dimension)* &
         (coupling*normalized_bracket_factor)**(total_length-1)
    associate(workspace => this%workspaces(total_length), &
              plan => this%recoupling%plans(total_length))
      current_buffer = 1
      call clear_terminal_layer(plan%max_inversions, current_buffer)
      call accumulate_terminal_seed_layer(plan%max_inversions, current_buffer)

      do inversion = plan%max_inversions, 1, -1
        next_buffer = 3-current_buffer
        call clear_terminal_layer(inversion-1, next_buffer)
        call accumulate_terminal_seed_layer(inversion-1, next_buffer)

        do layer_entry = plan%layer_node_start(inversion), &
                         plan%layer_node_start(inversion+1)-1
          node = plan%layer_nodes(layer_entry)
          swap_position = plan%nodes(node)%swap_position
          target_node = plan%nodes(node)%next_node
          if (swap_position < 1 .or. target_node < 1) &
               call fail('incomplete terminal recoupling ordering edge')
          source_slot = plan%nodes(node)%layer_slot
          target_slot = plan%nodes(target_node)%layer_slot
          call apply_sparse_scalar_add( &
               this%recoupling%swaps(total_length)%matrices(swap_position), &
               workspace%scalar(:, source_slot, current_buffer), &
               workspace%scalar(:, target_slot, next_buffer), &
               plan%nodes(node)%active_path, &
               workspace%initialized(:, target_slot, next_buffer), &
               plan%nodes(node)%destination_preinitialized)
        end do
        current_buffer = next_buffer
      end do

      canonical_slot = plan%nodes(plan%canonical_node)%layer_slot
      basis_amplitudes = (0.0_dp, 0.0_dp)
      do active_index = 1, &
           size(plan%nodes(plan%canonical_node)%active_path)
        path = plan%nodes(plan%canonical_node)%active_path(active_index)
        basis_amplitudes(path) = amplitude_scale* &
             workspace%scalar(path, canonical_slot, current_buffer)
      end do
    end associate

  contains

    subroutine clear_terminal_layer(inversion_layer, buffer)
      integer, intent(in) :: inversion_layer, buffer

      integer :: active, entry_layer, path_index, path_slot, workspace_node

      associate(workspace => this%workspaces(total_length), &
                plan => this%recoupling%plans(total_length))
        do entry_layer = plan%layer_node_start(inversion_layer), &
                         plan%layer_node_start(inversion_layer+1)-1
          workspace_node = plan%layer_nodes(entry_layer)
          path_slot = plan%nodes(workspace_node)%layer_slot
          if (size(plan%nodes(workspace_node)%active_path) == &
              plan%number_of_paths) then
            workspace%initialized(:, path_slot, buffer) = .false.
          else
            do active = 1, size(plan%nodes(workspace_node)%active_path)
              path_index = plan%nodes(workspace_node)%active_path(active)
              workspace%initialized(path_index, path_slot, buffer) = .false.
            end do
          end if
        end do
      end associate
    end subroutine clear_terminal_layer

    subroutine accumulate_terminal_seed_layer(inversion_layer, buffer)
      integer, intent(in) :: inversion_layer, buffer

      complex(dp) :: contracted, vertex_gluon(4), vertex_piece(4)
      integer :: initial_node, initial_slot, index_partition, layer_seed
      integer :: left_active, left_index, left_path, left_size, parent_path
      integer :: right_active, right_index, right_path
      integer :: signed_parent_path

      associate(workspace => this%workspaces(total_length), &
                evaluation => this%evaluations(pattern), &
                plan => this%recoupling%plans(total_length))
        do layer_seed = plan%layer_partition_start(inversion_layer), &
                        plan%layer_partition_start(inversion_layer+1)-1
          index_partition = plan%layer_partitions(layer_seed)
          left_index = int(evaluation%left_pattern(index_partition))
          right_index = int(evaluation%right_pattern(index_partition))
          left_size = plan%partitions(index_partition)%left_size
          initial_node = plan%partitions(index_partition)%initial_node
          initial_slot = plan%nodes(initial_node)%layer_slot
          associate(seed => plan%seeds(left_size), &
                    left => this%currents(left_index), &
                    right => this%currents(right_index), &
                    left_support => &
                         this%recoupling%canonical_support(left_size), &
                    right_support => &
                         this%recoupling%canonical_support( &
                              total_length-left_size))
            do right_active = 1, size(right_support%active_path)
              right_path = right_support%active_path(right_active)
              do left_active = 1, size(left_support%active_path)
                left_path = left_support%active_path(left_active)
                signed_parent_path = &
                     seed%signed_parent_path(left_path, right_path)
                parent_path = abs(signed_parent_path)

                call three_gluon(left%gluon(:, left_path), left%momentum, &
                                 right%gluon(:, right_path), right%momentum, &
                                 vertex_gluon)
                if (left%length >= 2) then
                  call aux_tensor_gluon_to_gluon( &
                       left%tensor(:, left_path), right%gluon(:, right_path), &
                       vertex_piece)
                  vertex_gluon = vertex_gluon+vertex_piece
                end if
                if (right%length >= 2) then
                  call gluon_aux_tensor_to_gluon( &
                       left%gluon(:, left_path), right%tensor(:, right_path), &
                       vertex_piece)
                  vertex_gluon = vertex_gluon+vertex_piece
                end if

                contracted = lorentz_contract(vertex_gluon, &
                                              closing_wavefunction)
                if (.not. workspace%initialized( &
                    parent_path, initial_slot, buffer)) then
                  if (signed_parent_path > 0) then
                    workspace%scalar(parent_path, initial_slot, buffer) = &
                         contracted
                  else
                    workspace%scalar(parent_path, initial_slot, buffer) = &
                         -contracted
                  end if
                  workspace%initialized(parent_path, initial_slot, buffer) = &
                       .true.
                else
                  if (signed_parent_path > 0) then
                    workspace%scalar(parent_path, initial_slot, buffer) = &
                         workspace%scalar(parent_path, initial_slot, buffer) + &
                         contracted
                  else
                    workspace%scalar(parent_path, initial_slot, buffer) = &
                         workspace%scalar(parent_path, initial_slot, buffer) - &
                         contracted
                  end if
                end if
              end do
            end do
          end associate
        end do
      end associate
    end subroutine accumulate_terminal_seed_layer

  end subroutine evaluate_terminal_subset

  subroutine validate_event(this, momenta, helicities)
    class(ampligluon_multiplet_t), intent(in) :: this
    real(dp), intent(in) :: momenta(0:, :)
    integer, intent(in) :: helicities(:)

    real(dp) :: mass_squared, momentum_scale, residual(0:3)
    integer :: external_leg

    if (size(momenta, 1) /= 4 .or. size(momenta, 2) /= this%total_gluons) &
         call fail('momentum array has the wrong shape')
    if (size(helicities) /= this%total_gluons) &
         call fail('helicity array has the wrong length')
    if (.not. all(ieee_is_finite(momenta))) call fail('momenta are not finite')
    if (any(abs(helicities) /= 1)) call fail('all gluon helicities must be -1 or +1')

    momentum_scale = max(1.0_dp, maxval(abs(momenta)))
    do external_leg = 1, this%total_gluons
      if (momenta(0, external_leg) <= 0.0_dp) &
           call fail('physical external energies must be positive')
      mass_squared = momenta(0, external_leg)**2 - &
           sum(momenta(1:3, external_leg)**2)
      if (abs(mass_squared) > momentum_tolerance*momentum_scale**2) &
           call fail('external momentum is not massless')
    end do
    residual = momenta(:, 1)+momenta(:, 2) - &
               sum(momenta(:, 3:this%total_gluons), dim=2)
    if (maxval(abs(residual)) > momentum_tolerance*momentum_scale) &
         call fail('external momenta do not conserve four-momentum')
  end subroutine validate_event

  integer function number_of_basis_amplitudes(this) result(number)
    class(ampligluon_multiplet_t), intent(in) :: this

    if (.not. this%initialized) call fail('amplitude object is not initialized')
    number = this%catalog%spaces(this%current_gluons)%number_of_paths
  end function number_of_basis_amplitudes

  integer function number_of_final_gluons(this) result(number)
    class(ampligluon_multiplet_t), intent(in) :: this

    if (.not. this%initialized) call fail('amplitude object is not initialized')
    number = this%total_gluons-2
  end function number_of_final_gluons

end module ampligluon_multiplet
