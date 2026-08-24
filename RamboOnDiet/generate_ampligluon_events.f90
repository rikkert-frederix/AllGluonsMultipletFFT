program generate_ampligluon_events
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: error_unit, int64
  use rambo_on_diet, only: generate_mom_from_x_v1
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: maximum_number_of_points = 999999
  integer, parameter :: maximum_total_gluons = 20
  real(dp), parameter :: default_centre_of_mass_energy = 1000.0_dp
  integer(int64), parameter :: default_seed = 1_int64
  real(dp), parameter :: random_boundary = 1.0e-12_dp

  character(len=4096) :: allocation_message, argument, event_filename
  character(len=4096) :: output_prefix
  integer, allocatable :: helicities(:, :)
  integer :: allocation_status, argument_count, final_gluons, ios
  integer :: number_of_points, point, total_gluons
  integer(int64) :: number_of_helicity_configurations, seed
  real(dp), allocatable :: masses(:), momenta(:, :), outgoing(:, :), x(:)
  real(dp) :: centre_of_mass_energy, phase_space_weight

  argument_count = command_argument_count()
  if (argument_count < 3 .or. argument_count > 5) call usage()

  call get_command_argument(1, argument)
  read(argument, *, iostat=ios) final_gluons
  if (ios /= 0 .or. final_gluons < 2) &
       call fail('FINAL_GLUONS must be an integer of at least two')

  call get_command_argument(2, argument)
  read(argument, *, iostat=ios) number_of_points
  if (ios /= 0 .or. number_of_points < 1 .or. &
       number_of_points > maximum_number_of_points) &
       call fail('POINTS must be between 1 and 999999')

  call get_command_argument(3, output_prefix)
  if (len_trim(output_prefix) == 0) call fail('OUTPUT_PREFIX must not be empty')
  if (len_trim(output_prefix)+len('_000000.event') > len(event_filename)) &
       call fail('OUTPUT_PREFIX is too long')

  centre_of_mass_energy = default_centre_of_mass_energy
  if (argument_count >= 4) then
    call get_command_argument(4, argument)
    read(argument, *, iostat=ios) centre_of_mass_energy
    if (ios /= 0 .or. .not. ieee_is_finite(centre_of_mass_energy) .or. &
         centre_of_mass_energy <= 0.0_dp) &
         call fail('SQRT_S must be finite and positive')
  end if

  seed = default_seed
  if (argument_count == 5) then
    call get_command_argument(5, argument)
    read(argument, *, iostat=ios) seed
    if (ios /= 0) call fail('SEED must be an integer')
  end if

  if (final_gluons > huge(final_gluons)-2) &
       call fail('FINAL_GLUONS is too large')
  total_gluons = final_gluons+2
  if (total_gluons > maximum_total_gluons) &
       call fail('exhaustive helicities support at most 20 total gluons')
  if (total_gluons >= bit_size(number_of_helicity_configurations)-1) &
       call fail('the number of helicity configurations overflows int64')
  number_of_helicity_configurations = shiftl(1_int64, total_gluons)
  if (number_of_helicity_configurations > int(huge(total_gluons), int64)) &
       call fail('the number of helicity configurations is too large')

  allocate(masses(final_gluons), x(3*final_gluons-4), &
           outgoing(0:3, final_gluons), momenta(0:3, total_gluons), &
           helicities(total_gluons, number_of_helicity_configurations), &
           stat=allocation_status, errmsg=allocation_message)
  if (allocation_status /= 0) &
       call fail('cannot allocate event data: '//trim(allocation_message))
  masses = 0.0_dp
  call fill_helicities(helicities)
  call initialize_random_seed(seed)

  momenta = 0.0_dp
  momenta(:, 1) = [centre_of_mass_energy/2.0_dp, 0.0_dp, 0.0_dp, &
                   centre_of_mass_energy/2.0_dp]
  momenta(:, 2) = [centre_of_mass_energy/2.0_dp, 0.0_dp, 0.0_dp, &
                  -centre_of_mass_energy/2.0_dp]

  write(*, '(a,1x,i0)') 'SEED', seed
  do point = 1, number_of_points
    call random_number(x)
    x = max(random_boundary, min(1.0_dp-random_boundary, x))
    call generate_mom_from_x_v1(final_gluons, x, masses, &
         centre_of_mass_energy, outgoing, phase_space_weight)
    momenta(:, 3:total_gluons) = outgoing
    call validate_phase_space(momenta, phase_space_weight)

    write(argument, '(i6.6)') point
    event_filename = trim(output_prefix)//'_'//trim(argument)//'.event'
    call write_event(trim(event_filename), final_gluons, momenta, helicities)
    write(*, '(a,1x,a,1x,a,1x,es25.17e3)') &
         'EVENT', trim(event_filename), 'WEIGHT', phase_space_weight
  end do

contains

  subroutine usage()
    write(error_unit, '(a)') &
         'usage: generate_ampligluon_events '// &
         'FINAL_GLUONS POINTS OUTPUT_PREFIX [SQRT_S [SEED]]'
    error stop 2
  end subroutine usage

  subroutine fail(message)
    character(len=*), intent(in) :: message

    write(error_unit, '(a)') 'generate_ampligluon_events: '//trim(message)
    error stop 1
  end subroutine fail

  subroutine initialize_random_seed(user_seed)
    integer(int64), intent(in) :: user_seed

    integer, allocatable :: seed_values(:)
    integer :: seed_index, seed_size
    integer(int64) :: state
    integer(int64), parameter :: multiplier = 48271_int64
    integer(int64), parameter :: modulus = 2147483647_int64

    call random_seed(size=seed_size)
    allocate(seed_values(seed_size))
    state = modulo(user_seed, modulus-1_int64)+1_int64
    do seed_index = 1, seed_size
      state = modulo(multiplier*state, modulus)
      seed_values(seed_index) = int(state)
    end do
    call random_seed(put=seed_values)
  end subroutine initialize_random_seed

  subroutine fill_helicities(configurations)
    integer, intent(out) :: configurations(:, :)

    integer :: leg
    integer(int64) :: configuration

    do configuration = 0_int64, size(configurations, 2, kind=int64)-1_int64
      do leg = 1, size(configurations, 1)
        if (btest(configuration, size(configurations, 1)-leg)) then
          configurations(leg, configuration+1_int64) = 1
        else
          configurations(leg, configuration+1_int64) = -1
        end if
      end do
    end do
  end subroutine fill_helicities

  subroutine validate_phase_space(event_momenta, weight)
    real(dp), intent(in) :: event_momenta(0:, :), weight

    real(dp) :: mass_squared, momentum_scale, residual(0:3)
    integer :: leg

    if (.not. all(ieee_is_finite(event_momenta))) &
         call fail('RamboOnDiet produced non-finite momenta')
    if (.not. ieee_is_finite(weight) .or. weight <= 0.0_dp) &
         call fail('RamboOnDiet produced an invalid phase-space weight')

    momentum_scale = max(1.0_dp, maxval(abs(event_momenta)))
    do leg = 1, size(event_momenta, 2)
      if (event_momenta(0, leg) <= 0.0_dp) &
           call fail('RamboOnDiet produced a non-positive energy')
      mass_squared = event_momenta(0, leg)**2- &
           sum(event_momenta(1:3, leg)**2)
      if (abs(mass_squared) > 1.0e-10_dp*momentum_scale**2) &
           call fail('RamboOnDiet produced a non-massless momentum')
    end do

    residual = event_momenta(:, 1)+event_momenta(:, 2)- &
         sum(event_momenta(:, 3:size(event_momenta, 2)), dim=2)
    if (maxval(abs(residual)) > 1.0e-10_dp*momentum_scale) &
         call fail('RamboOnDiet momenta do not conserve four-momentum')
  end subroutine validate_phase_space

  subroutine write_event(filename, number_of_final_gluons, &
                         event_momenta, configurations)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: number_of_final_gluons
    real(dp), intent(in) :: event_momenta(0:, :)
    integer, intent(in) :: configurations(:, :)

    integer :: configuration, leg, output_unit, status

    open(newunit=output_unit, file=filename, status='replace', action='write', &
         form='formatted', iostat=status)
    if (status /= 0) call fail('cannot create event file: '//trim(filename))

    write(output_unit, '(a)') 'AMPLIGLUON_EVENT_V1'
    write(output_unit, '(a,1x,i0)') 'FINAL_GLUONS', number_of_final_gluons
    write(output_unit, '(a,1x,es25.17e3)') 'STRONG_COUPLING', 1.0_dp
    write(output_unit, '(a)') 'BEGIN_MOMENTA'
    do leg = 1, size(event_momenta, 2)
      write(output_unit, '(*(es25.17e3,1x))') event_momenta(:, leg)
    end do
    write(output_unit, '(a)') 'END_MOMENTA'
    write(output_unit, '(a,1x,i0)') 'NHELICITIES', size(configurations, 2)
    write(output_unit, '(a)') 'BEGIN_HELICITIES'
    do configuration = 1, size(configurations, 2)
      write(output_unit, '(*(sp,i0,1x))') configurations(:, configuration)
    end do
    write(output_unit, '(a)') 'END_HELICITIES'
    close(output_unit, iostat=status)
    if (status /= 0) call fail('cannot close event file: '//trim(filename))
  end subroutine write_event

end program generate_ampligluon_events
