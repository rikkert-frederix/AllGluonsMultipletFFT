module multiplet_mhv
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use ampligluon_multiplet_kinds, only: dp, fail
  use multiplet_paths, only: path_catalog_t
  use multiplet_radiation, only: multiplet_radiation_system_t, &
       multiplet_radiation_workspace_t
  use qcd_kinematics, only: aux_tensor_gluon_to_gluon, &
       external_massless_vector, gluon_aux_tensor_to_gluon, &
       lorentz_contract, massless_vector_propagator, three_gluon, &
       two_gluon_to_aux_tensor
  use recoupling_plan, only: recoupling_system_t, sparse_matrix_t
  implicit none
  private

  real(dp), parameter :: adjoint_dimension = 8.0_dp
  real(dp), parameter :: normalized_bracket_factor = sqrt(6.0_dp)
  real(dp), parameter :: spinor_safety_factor = 4096.0_dp
  real(dp), parameter :: minimum_pair_separation = 2.0e-3_dp

  ! MHV and anti-MHV amplitudes are built directly in the sequential,
  ! orthonormal multiplet path basis.  The colour recursion is an inverse-soft
  ! sum of normalized radiation maps.  No ordered-colour vector or colour Gram
  ! matrix is constructed.  One ordered Berends--Giele amplitude fixes only
  ! the momentum-chart phase of the analytic spinors relative to the public
  ! polarization convention.
  type, public :: multiplet_mhv_t
    private
    logical :: initialized = .false.
    integer :: total_gluons = 0
    integer :: current_gluons = 0
    integer :: antisymmetric_seed_path = 0
    type(multiplet_radiation_system_t) :: radiation
    type(multiplet_radiation_workspace_t) :: radiation_workspace
    complex(dp), allocatable :: angle(:, :)
    complex(dp), allocatable :: square(:, :)
    complex(dp), allocatable :: current(:)
    complex(dp), allocatable :: next(:)
    complex(dp), allocatable :: compensation(:)
    complex(dp), allocatable :: contribution(:)
    complex(dp), allocatable :: bg_gluon(:, :, :)
    complex(dp), allocatable :: bg_tensor(:, :, :)
    real(dp), allocatable :: bg_momentum(:, :, :)
    integer, allocatable :: active_order(:)
  contains
    procedure, public :: evaluate => evaluate_multiplet_mhv
    procedure, public :: initialize => initialize_multiplet_mhv
  end type multiplet_mhv_t

contains

  subroutine initialize_multiplet_mhv(this, catalog)
    class(multiplet_mhv_t), intent(inout) :: this
    type(path_catalog_t), intent(in) :: catalog

    integer :: maximum_paths, path

    if (.not. allocated(catalog%spaces) .or. catalog%max_length < 3) &
         call fail('MHV evaluator needs multiplet paths through length three')
    this%initialized = .false.
    this%current_gluons = catalog%max_length
    this%total_gluons = this%current_gluons+1
    maximum_paths = catalog%spaces(this%current_gluons)%number_of_paths

    this%antisymmetric_seed_path = 0
    do path = 1, catalog%spaces(2)%number_of_paths
      if (catalog%spaces(2)%rep_p(2, path) /= 1 .or. &
          catalog%spaces(2)%rep_q(2, path) /= 1 .or. &
          catalog%spaces(2)%multiplicity(2, path) /= 1) cycle
      if (this%antisymmetric_seed_path /= 0) &
           call fail('MHV seed has more than one antisymmetric adjoint path')
      this%antisymmetric_seed_path = path
    end do
    if (this%antisymmetric_seed_path == 0) &
         call fail('MHV seed is missing the antisymmetric adjoint path')

    call this%radiation%build(catalog)
    if (allocated(this%angle)) deallocate(this%angle, this%square)
    if (allocated(this%current)) then
      deallocate(this%current, this%next, this%compensation, this%contribution)
    end if
    if (allocated(this%bg_gluon)) then
      deallocate(this%bg_gluon, this%bg_tensor, this%bg_momentum)
    end if
    if (allocated(this%active_order)) deallocate(this%active_order)
    allocate(this%angle(this%total_gluons, this%total_gluons), &
             this%square(this%total_gluons, this%total_gluons))
    allocate(this%current(maximum_paths), this%next(maximum_paths), &
             this%compensation(maximum_paths), &
             this%contribution(maximum_paths))
    allocate(this%bg_gluon(4, this%current_gluons, this%current_gluons), &
             this%bg_tensor(6, this%current_gluons, this%current_gluons), &
             this%bg_momentum(0:3, this%current_gluons, this%current_gluons))
    allocate(this%active_order(this%total_gluons))
    this%initialized = .true.
  end subroutine initialize_multiplet_mhv

  subroutine evaluate_multiplet_mhv(this, recoupling, outgoing_momenta, &
                                    physical_helicities, coupling, &
                                    basis_amplitudes, succeeded)
    class(multiplet_mhv_t), intent(inout) :: this
    type(recoupling_system_t), intent(in) :: recoupling
    real(dp), intent(in) :: outgoing_momenta(0:, :)
    integer, intent(in) :: physical_helicities(:)
    real(dp), intent(in) :: coupling
    complex(dp), intent(out) :: basis_amplitudes(:)
    logical, intent(out) :: succeeded

    complex(dp) :: phase, scale, soft_weight
    complex(dp), allocatable :: bracket(:, :)
    integer :: emitted, emitter, negative_count, old_length
    integer :: positive_count, reference, source_size, special(2), target_size
    logical :: use_angle

    succeeded = .false.
    if (.not. this%initialized) call fail('MHV evaluator is not initialized')
    if (size(outgoing_momenta, 1) /= 4 .or. &
        size(outgoing_momenta, 2) /= this%total_gluons .or. &
        size(physical_helicities) /= this%total_gluons) &
         call fail('wrong event shape in MHV evaluator')
    if (size(basis_amplitudes) /= &
        this%radiation%target_dimension(this%current_gluons-1)) &
         call fail('wrong output basis size in MHV evaluator')

    negative_count = count(physical_helicities(1:2) > 0) + &
                     count(physical_helicities(3:) < 0)
    positive_count = this%total_gluons-negative_count
    if (negative_count == 2) then
      use_angle = .true.
      call find_special_legs(physical_helicities, -1, special)
    else if (positive_count == 2) then
      use_angle = .false.
      call find_special_legs(physical_helicities, 1, special)
    else
      return
    end if

    ! The inverse-soft representation is exact, but cancellations between
    ! its terms become ill-conditioned when any two external directions are
    ! extremely close.  The planar phase check sees only adjacent pairs and
    ! cannot diagnose every such channel. Use the dimensionless event-frame
    ! pair separation |<ij>|/sqrt(|E_i E_j|) (equal to its square-bracket
    ! counterpart for real momenta) and let the general multiplet BG recursion
    ! handle the deliberately conservative near-collinear region.
    if (.not. pair_separations_are_safe(outgoing_momenta)) return

    call massless_spinor_products(outgoing_momenta, this%angle, this%square, &
                                  succeeded)
    if (.not. succeeded) return
    allocate(bracket(this%total_gluons, this%total_gluons))
    if (use_angle) then
      bracket = this%angle
    else
      bracket = this%square
    end if

    reference = first_regular_leg(this%total_gluons, special)
    call build_active_order(this%total_gluons, special, reference, &
                            this%active_order)
    if (.not. inverse_soft_denominators_are_safe( &
        bracket, this%active_order, special(1))) then
      succeeded = .false.
      return
    end if

    this%current = (0.0_dp, 0.0_dp)
    this%current(this%antisymmetric_seed_path) = (1.0_dp, 0.0_dp)
    do emitted = 4, this%total_gluons
      old_length = emitted-2
      source_size = this%radiation%source_dimension(old_length)
      target_size = this%radiation%target_dimension(old_length)
      this%next(1:target_size) = (0.0_dp, 0.0_dp)
      this%compensation(1:target_size) = (0.0_dp, 0.0_dp)
      do emitter = 1, emitted-1
        if (this%active_order(emitter) == special(1)) cycle
        soft_weight = &
             bracket(this%active_order(emitter), special(1))/ &
             (bracket(this%active_order(emitter), this%active_order(emitted))* &
              bracket(this%active_order(emitted), special(1)))
        call this%radiation%apply(recoupling, old_length, emitter, &
             this%current(1:source_size), this%contribution(1:target_size), &
             this%radiation_workspace)
        call compensated_vector_add(soft_weight, &
             this%contribution(1:target_size), this%next(1:target_size), &
             this%compensation(1:target_size))
      end do
      this%current(1:target_size) = this%next(1:target_size)
    end do

    call reorder_to_canonical(this, recoupling, this%current(1:size(basis_amplitudes)))
    call ordered_phase(this, outgoing_momenta, physical_helicities, bracket, &
                       special, phase, succeeded)
    if (.not. succeeded) return
    scale = sqrt(adjoint_dimension)* &
         (coupling*normalized_bracket_factor)**(this%total_gluons-2)* &
         three_point_amplitude(bracket, special, reference)*phase
    basis_amplitudes = scale*this%current(1:size(basis_amplitudes))
    succeeded = all(ieee_is_finite(real(basis_amplitudes, dp))) .and. &
                all(ieee_is_finite(aimag(basis_amplitudes)))
  end subroutine evaluate_multiplet_mhv

  subroutine find_special_legs(physical_helicities, outgoing_sign, special)
    integer, intent(in) :: physical_helicities(:), outgoing_sign
    integer, intent(out) :: special(2)

    integer :: found, leg, outgoing_helicity

    found = 0
    do leg = 1, size(physical_helicities)
      outgoing_helicity = physical_helicities(leg)
      if (leg <= 2) outgoing_helicity = -outgoing_helicity
      if (outgoing_helicity /= outgoing_sign) cycle
      found = found+1
      if (found <= 2) special(found) = leg
    end do
    if (found /= 2) call fail('MHV special-leg count is inconsistent')
  end subroutine find_special_legs

  integer function first_regular_leg(number, special) result(leg)
    integer, intent(in) :: number, special(2)

    do leg = 1, number
      if (all(leg /= special)) return
    end do
    call fail('MHV seed has no regular leg')
  end function first_regular_leg

  subroutine build_active_order(number, special, reference, order)
    integer, intent(in) :: number, special(2), reference
    integer, intent(out) :: order(number)

    integer :: leg, position

    order(1:3) = [special, reference]
    position = 3
    do leg = 1, number
      if (any(leg == order(1:3))) cycle
      position = position+1
      order(position) = leg
    end do
    if (position /= number) call fail('incomplete MHV inverse-soft order')
  end subroutine build_active_order

  logical function pair_separations_are_safe(momenta) result(safe)
    real(dp), intent(in) :: momenta(0:, :)

    real(dp) :: denominator, dot_product, separation
    integer :: first, second

    safe = .true.
    do second = 2, size(momenta, 2)
      do first = 1, second-1
        denominator = abs(momenta(0, first)*momenta(0, second))
        if (denominator <= tiny(1.0_dp)) then
          safe = .false.
          return
        end if
        dot_product = momenta(0, first)*momenta(0, second) - &
             sum(momenta(1:3, first)*momenta(1:3, second))
        separation = sqrt(abs(2.0_dp*dot_product)/denominator)
        if (.not. ieee_is_finite(separation) .or. &
            separation < minimum_pair_separation) then
          safe = .false.
          return
        end if
      end do
    end do
  end function pair_separations_are_safe

  logical function inverse_soft_denominators_are_safe(bracket, order, reference) &
       result(safe)
    complex(dp), intent(in) :: bracket(:, :)
    integer, intent(in) :: order(:), reference

    real(dp) :: threshold
    integer :: emitted, emitter

    threshold = spinor_safety_factor*epsilon(1.0_dp)* &
                max(1.0_dp, maxval(abs(bracket)))
    safe = .true.
    do emitted = 4, size(order)
      if (abs(bracket(order(emitted), reference)) <= threshold) then
        safe = .false.
        return
      end if
      do emitter = 1, emitted-1
        if (order(emitter) == reference) cycle
        if (abs(bracket(order(emitter), order(emitted))) <= threshold) then
          safe = .false.
          return
        end if
      end do
    end do
  end function inverse_soft_denominators_are_safe

  subroutine compensated_vector_add(weight, source, destination, compensation)
    complex(dp), intent(in) :: weight, source(:)
    complex(dp), intent(inout) :: destination(:), compensation(:)

    complex(dp) :: corrected, updated
    integer :: entry

    if (size(destination) /= size(source) .or. &
        size(compensation) /= size(source)) &
         call fail('wrong vector shape in compensated MHV sum')
    do entry = 1, size(source)
      corrected = weight*source(entry)-compensation(entry)
      updated = destination(entry)+corrected
      compensation(entry) = (updated-destination(entry))-corrected
      destination(entry) = updated
    end do
  end subroutine compensated_vector_add

  subroutine reorder_to_canonical(this, recoupling, vector)
    class(multiplet_mhv_t), intent(inout) :: this
    type(recoupling_system_t), intent(in) :: recoupling
    complex(dp), intent(inout) :: vector(:)

    integer :: left, position, temporary

    do
      position = 0
      do left = 1, this%total_gluons-1
        if (this%active_order(left) <= this%active_order(left+1)) cycle
        position = left
        exit
      end do
      if (position == 0) exit
      call apply_swap(recoupling%swaps(this%current_gluons)%matrices(position), &
                      vector, this%next(1:size(vector)))
      vector = this%next(1:size(vector))
      temporary = this%active_order(position)
      this%active_order(position) = this%active_order(position+1)
      this%active_order(position+1) = temporary
    end do
  end subroutine reorder_to_canonical

  subroutine apply_swap(matrix, source, destination)
    type(sparse_matrix_t), intent(in) :: matrix
    complex(dp), intent(in) :: source(:)
    complex(dp), intent(out) :: destination(:)

    integer :: column, entry, signed_row

    if (size(source) /= matrix%size .or. size(destination) /= matrix%size) &
         call fail('wrong vector size in MHV basis swap')
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

  subroutine massless_spinor_products(momenta, angle, square, succeeded)
    real(dp), intent(in) :: momenta(0:, :)
    complex(dp), intent(out), target :: angle(:, :), square(:, :)
    logical, intent(out) :: succeeded

    complex(dp), allocatable :: lambda(:, :), lambda_tilde(:, :)
    complex(dp) :: root
    integer :: first, leg, second

    allocate(lambda(2, size(momenta, 2)), &
             lambda_tilde(2, size(momenta, 2)))
    succeeded = .false.
    do leg = 1, size(momenta, 2)
      if (abs(momenta(0, leg)+momenta(3, leg)) >= &
          abs(momenta(0, leg)-momenta(3, leg))) then
        root = sqrt(cmplx(momenta(0, leg)+momenta(3, leg), 0.0_dp, dp))
        if (abs(root) <= tiny(1.0_dp)) return
        lambda(1, leg) = root
        lambda(2, leg) = &
             cmplx(momenta(1, leg), momenta(2, leg), dp)/root
        lambda_tilde(1, leg) = root
        lambda_tilde(2, leg) = &
             cmplx(momenta(1, leg), -momenta(2, leg), dp)/root
      else
        root = sqrt(cmplx(momenta(0, leg)-momenta(3, leg), 0.0_dp, dp))
        if (abs(root) <= tiny(1.0_dp)) return
        lambda(1, leg) = &
             cmplx(momenta(1, leg), -momenta(2, leg), dp)/root
        lambda(2, leg) = root
        lambda_tilde(1, leg) = &
             cmplx(momenta(1, leg), momenta(2, leg), dp)/root
        lambda_tilde(2, leg) = root
      end if
    end do
    do second = 1, size(momenta, 2)
      do first = 1, size(momenta, 2)
        angle(first, second) = &
             lambda(1, first)*lambda(2, second)- &
             lambda(2, first)*lambda(1, second)
        square(first, second) = &
             lambda_tilde(1, first)*lambda_tilde(2, second)- &
             lambda_tilde(2, first)*lambda_tilde(1, second)
      end do
    end do
    succeeded = all(ieee_is_finite(real(angle, dp))) .and. &
                all(ieee_is_finite(aimag(angle))) .and. &
                all(ieee_is_finite(real(square, dp))) .and. &
                all(ieee_is_finite(aimag(square)))
  end subroutine massless_spinor_products

  pure complex(dp) function three_point_amplitude(bracket, special, reference) &
       result(amplitude)
    complex(dp), intent(in) :: bracket(:, :)
    integer, intent(in) :: special(2), reference

    amplitude = bracket(special(1), special(2))**3/ &
         (bracket(special(2), reference)*bracket(reference, special(1)))
  end function three_point_amplitude

  subroutine ordered_phase(this, outgoing_momenta, physical_helicities, &
                           bracket, special, phase, succeeded)
    class(multiplet_mhv_t), intent(inout) :: this
    real(dp), intent(in) :: outgoing_momenta(0:, :)
    integer, intent(in) :: physical_helicities(:), special(2)
    complex(dp), intent(in) :: bracket(:, :)
    complex(dp), intent(out) :: phase
    logical, intent(out) :: succeeded

    complex(dp) :: ordered_amplitude, parke_taylor, product
    real(dp) :: magnitude, threshold
    integer :: leg

    product = (1.0_dp, 0.0_dp)
    do leg = 1, this%total_gluons-1
      product = product*bracket(leg, leg+1)
    end do
    product = product*bracket(this%total_gluons, 1)
    threshold = spinor_safety_factor*epsilon(1.0_dp)* &
                max(1.0_dp, maxval(abs(bracket)))**this%total_gluons
    if (abs(product) <= threshold) then
      succeeded = .false.
      return
    end if
    parke_taylor = bracket(special(1), special(2))**4/product
    call evaluate_one_order(this, outgoing_momenta, physical_helicities, &
                            ordered_amplitude)
    if (abs(parke_taylor) <= tiny(1.0_dp)) then
      succeeded = .false.
      return
    end if
    phase = ordered_amplitude/parke_taylor
    magnitude = abs(phase)
    if (.not. ieee_is_finite(magnitude) .or. magnitude <= tiny(1.0_dp) .or. &
        abs(magnitude-1.0_dp) > 2.0e-8_dp) then
      succeeded = .false.
      return
    end if
    succeeded = .true.
  end subroutine ordered_phase

  subroutine evaluate_one_order(this, outgoing_momenta, physical_helicities, &
                                amplitude)
    class(multiplet_mhv_t), intent(inout) :: this
    real(dp), intent(in) :: outgoing_momenta(0:, :)
    integer, intent(in) :: physical_helicities(:)
    complex(dp), intent(out) :: amplitude

    complex(dp) :: closing_wavefunction(4), tensor_vertex(6), vector_vertex(4)
    integer :: first, last, left_last, leg, length

    this%bg_gluon = (0.0_dp, 0.0_dp)
    this%bg_tensor = (0.0_dp, 0.0_dp)
    this%bg_momentum = 0.0_dp
    do leg = 1, this%current_gluons
      this%bg_momentum(:, leg, leg) = outgoing_momenta(:, leg)
      call external_massless_vector(outgoing_momenta(:, leg), &
           physical_helicities(leg), this%bg_gluon(:, leg, leg))
    end do
    call external_massless_vector(outgoing_momenta(:, this%total_gluons), &
         physical_helicities(this%total_gluons), closing_wavefunction)

    do length = 2, this%current_gluons
      do first = 1, this%current_gluons-length+1
        last = first+length-1
        this%bg_momentum(:, first, last) = &
             this%bg_momentum(:, first, last-1)+ &
             this%bg_momentum(:, last, last)
        do left_last = first, last-1
          call three_gluon(this%bg_gluon(:, first, left_last), &
               this%bg_momentum(:, first, left_last), &
               this%bg_gluon(:, left_last+1, last), &
               this%bg_momentum(:, left_last+1, last), vector_vertex)
          this%bg_gluon(:, first, last) = &
               this%bg_gluon(:, first, last)+vector_vertex
          if (left_last > first) then
            call aux_tensor_gluon_to_gluon( &
                 this%bg_tensor(:, first, left_last), &
                 this%bg_gluon(:, left_last+1, last), vector_vertex)
            this%bg_gluon(:, first, last) = &
                 this%bg_gluon(:, first, last)+vector_vertex
          end if
          if (left_last+1 < last) then
            call gluon_aux_tensor_to_gluon( &
                 this%bg_gluon(:, first, left_last), &
                 this%bg_tensor(:, left_last+1, last), vector_vertex)
            this%bg_gluon(:, first, last) = &
                 this%bg_gluon(:, first, last)+vector_vertex
          end if
          if (length < this%current_gluons) then
            call two_gluon_to_aux_tensor( &
                 this%bg_gluon(:, first, left_last), &
                 this%bg_gluon(:, left_last+1, last), tensor_vertex)
            this%bg_tensor(:, first, last) = &
                 this%bg_tensor(:, first, last)+tensor_vertex
          end if
        end do
        if (length < this%current_gluons) call massless_vector_propagator( &
             this%bg_gluon(:, first:first, last), &
             this%bg_momentum(:, first, last))
      end do
    end do
    amplitude = lorentz_contract( &
         this%bg_gluon(:, 1, this%current_gluons), closing_wavefunction)
  end subroutine evaluate_one_order

end module multiplet_mhv
