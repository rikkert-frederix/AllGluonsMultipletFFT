program test_ampligluon_adjoint
  use adjoint_colour_matrix, only: adjoint_colour_matrix_t
  use ampligluon_adjoint, only: ampligluon_adjoint_t
  use ampligluon_common, only: dp, fail
  implicit none

  call check_multiplicity(4, 2, 3456.0_dp, 216.0_dp)
  call check_multiplicity(5, 6, 0.839808_dp, 0.010368_dp)
  call check_multiplicity(6, 24, 9.3333333333333397e-4_dp, &
       1.2049861549399821e-5_dp)
  call check_multiplicity(7, 120, 6.9826518403237469e-6_dp, &
       5.8326287358333135e-9_dp)
  call check_multiplicity(8, 720, 2.0200927180151527e-8_dp, &
       5.8732224850953953e-10_dp)
  call check_multiplicity(9, 5040, 2.9947601653227329e-8_dp, &
       1.2609237070482274e-11_dp)
  call check_multiplicity(10, 40320, 5.0926380699009548e-10_dp, &
       2.3105063855067617e-13_dp)
  call check_mhv_switch()
  call check_colour_matrix()
  call check_colour_contraction_modes()
  call check_public_colour_contraction_modes()
  call check_options()
  write(*, '(a)') 'AmpliGluonAdjoint regression: PASS'

contains

  subroutine check_colour_matrix()
    type(adjoint_colour_matrix_t) :: colour
    complex(dp), allocatable :: vector(:)
    integer, allocatable :: permutations(:, :)
    integer :: column
    integer, parameter :: first_row(6) = [1728, 864, 864, 432, 432, 0]
    real(dp) :: result

    call colour%initialize(4, permutations)
    allocate(vector(2))
    vector = (0.0_dp, 0.0_dp)
    vector(1) = (1.0_dp, 0.0_dp)
    call contract_readonly(colour, vector, result)
    call require_close(result, 288.0_dp, 'four-gluon Gram diagonal')
    vector(2) = (1.0_dp, 0.0_dp)
    call colour%contract(vector, result)
    call require_close(result, 864.0_dp, 'four-gluon Gram off-diagonal')

    call colour%initialize(5, permutations)
    deallocate(vector)
    allocate(vector(6))
    do column = 1, 6
      vector = (0.0_dp, 0.0_dp)
      vector(1) = (1.0_dp, 0.0_dp)
      if (column == 1) then
        call colour%contract(vector, result)
      else
        vector(column) = (1.0_dp, 0.0_dp)
        call colour%contract(vector, result)
        result = 0.5_dp*(result-3456.0_dp)
      end if
      call require_close(result, real(first_row(column), dp), &
           'five-gluon Gram first row')
    end do

    ! This high-multiplicity colour-only check guards against restoring the
    ! former quadratic packed Gram matrix, which could not initialize at N=10.
    call colour%initialize(10, permutations)
    deallocate(vector)
    allocate(vector(40320))
    vector = (0.0_dp, 0.0_dp)
    vector(1) = (1.0_dp, 0.0_dp)
    call colour%contract(vector, result)
    call require_close(result, 13436928.0_dp, 'ten-gluon Gram diagonal')
  end subroutine check_colour_matrix

  subroutine check_colour_contraction_modes()
    type(adjoint_colour_matrix_t) :: direct_colour, fourier_colour
    complex(dp), allocatable :: vector(:)
    integer, allocatable :: permutations(:, :)
    integer :: order, row, total_gluons
    real(dp) :: direct, direct_reused, fourier, scale

    do total_gluons = 4, 7
      call fourier_colour%initialize(total_gluons, permutations, &
           use_colour_fft=.true.)
      order = fourier_colour%number_of_basis_elements()
      call require(fourier_colour%number_of_stored_kernel_values() == order, &
           'Fourier adjoint kernel storage is not linear')
      call direct_colour%initialize(total_gluons, permutations, &
           use_colour_fft=.false.)
      call require(direct_colour%number_of_basis_elements() == order, &
           'direct and Fourier adjoint dimensions disagree')
      call require(direct_colour%number_of_stored_kernel_values() == order, &
           'direct adjoint kernel storage is not linear')
      if (allocated(vector)) deallocate(vector)
      allocate(vector(order))

      do row = 1, order
        vector(row) = cmplx( &
             sin(0.173_dp*real(row, dp))+0.003_dp*real(row, dp), &
             cos(0.119_dp*real(row, dp))-0.002_dp*real(row, dp), dp)
      end do
      call fourier_colour%contract(vector, fourier)
      call direct_colour%contract(vector, direct)
      call direct_colour%contract_reusing_workspace(vector, direct_reused)
      scale = max(1.0_dp, abs(fourier))
      call require(abs(direct-fourier) <= 2.0e-10_dp*scale, &
           'direct and Fourier adjoint contractions disagree')
      call require(abs(direct_reused-direct) <= 2.0e-12_dp*scale, &
           'direct reusable adjoint contraction disagrees')

      do row = 1, order
        vector(row) = cmplx( &
             cos(0.151_dp*real(row, dp))-0.004_dp*real(row, dp), &
             sin(0.107_dp*real(row, dp))+0.006_dp*real(row, dp), dp)
      end do
      call fourier_colour%contract_reusing_workspace(vector, fourier)
      call direct_colour%contract(vector, direct)
      scale = max(1.0_dp, abs(fourier))
      call require(abs(direct-fourier) <= 2.0e-10_dp*scale, &
           'reused direct and Fourier adjoint contractions disagree')
    end do

    ! Reinitialize one object across both representations.  This exercises
    ! release and replacement of the representation-specific storage.
    call direct_colour%initialize(6, permutations, use_colour_fft=.true.)
    call direct_colour%initialize(6, permutations, use_colour_fft=.false.)
    call direct_colour%contract(vector(1:24), direct)
    call fourier_colour%initialize(6, permutations, use_colour_fft=.true.)
    call fourier_colour%contract(vector(1:24), fourier)
    call require_close(direct, fourier, &
         'adjoint colour mode reinitialization mismatch')
  end subroutine check_colour_contraction_modes

  subroutine check_public_colour_contraction_modes()
    type(ampligluon_adjoint_t) :: direct_amplitude, fourier_amplitude
    real(dp) :: direct, fourier, momenta(0:3, 6)
    integer :: helicities(6)

    call fill_momenta(6, momenta)
    call direct_amplitude%initialize(4, use_colour_fft=.false.)
    call fourier_amplitude%initialize(4)

    helicities = -1
    call compare_public_colour_modes(direct_amplitude, fourier_amplitude, &
         momenta, helicities, .true., 'MHV colour mode mismatch')
    call compare_public_colour_modes(direct_amplitude, fourier_amplitude, &
         momenta, helicities, .false., &
         'forced-BG colour mode mismatch')

    helicities = [-1, 1, -1, 1, -1, 1]
    call compare_public_colour_modes(direct_amplitude, fourier_amplitude, &
         momenta, helicities, .true., &
         'non-MHV colour mode mismatch')

    ! The choice may be changed by reinitializing an existing amplitude.
    call direct_amplitude%initialize(4, use_colour_fft=.true.)
    call direct_amplitude%evaluate(momenta, helicities, direct, &
         strong_coupling=0.73_dp, average_initial_colours=.true.)
    call fourier_amplitude%evaluate(momenta, helicities, fourier, &
         strong_coupling=0.73_dp, average_initial_colours=.true.)
    call require_close(direct, fourier, &
         'public adjoint colour mode reinitialization mismatch')

  end subroutine check_public_colour_contraction_modes

  subroutine compare_public_colour_modes(direct_amplitude, fourier_amplitude, &
       momenta, helicities, use_analytic_mhv, description)
    type(ampligluon_adjoint_t), intent(inout) :: direct_amplitude
    type(ampligluon_adjoint_t), intent(inout) :: fourier_amplitude
    real(dp), intent(in) :: momenta(0:, :)
    integer, intent(in) :: helicities(:)
    logical, intent(in) :: use_analytic_mhv
    character(len=*), intent(in) :: description

    complex(dp), allocatable :: direct_basis(:), fourier_basis(:)
    real(dp) :: direct, fourier

    call direct_amplitude%evaluate(momenta, helicities, direct, &
         strong_coupling=0.73_dp, average_initial_colours=.true., &
         adjoint_amplitudes=direct_basis, &
         use_analytic_mhv=use_analytic_mhv)
    call fourier_amplitude%evaluate(momenta, helicities, fourier, &
         strong_coupling=0.73_dp, average_initial_colours=.true., &
         adjoint_amplitudes=fourier_basis, &
         use_analytic_mhv=use_analytic_mhv)
    call require_close(direct, fourier, description)
    call require_complex_vector_close(direct_basis, fourier_basis, &
         trim(description)//' basis amplitudes')
  end subroutine compare_public_colour_modes

  subroutine contract_readonly(colour, vector, result)
    type(adjoint_colour_matrix_t), intent(in) :: colour
    complex(dp), intent(in) :: vector(:)
    real(dp), intent(out) :: result

    call colour%contract(vector, result)
  end subroutine contract_readonly

  subroutine check_multiplicity(total_gluons, expected_dimension, &
                                expected_all_minus, expected_alternating)
    integer, intent(in) :: total_gluons, expected_dimension
    real(dp), intent(in) :: expected_all_minus, expected_alternating

    type(ampligluon_adjoint_t) :: amplitude
    complex(dp), allocatable :: analytic_basis(:), basis(:), recursive_basis(:)
    integer, allocatable :: helicities(:)
    real(dp), allocatable :: momenta(:, :)
    real(dp) :: analytic_matrix2, matrix2, recursive_matrix2
    integer :: leg

    allocate(helicities(total_gluons), momenta(0:3, total_gluons))
    call fill_momenta(total_gluons, momenta)
    call amplitude%initialize(total_gluons-2)
    call require(amplitude%number_of_final_gluons() == total_gluons-2, &
         'wrong final-state count')
    call require(amplitude%number_of_basis_amplitudes() == expected_dimension, &
         'wrong adjoint-basis dimension')

    helicities = -1
    call amplitude%evaluate(momenta, helicities, matrix2, &
         adjoint_amplitudes=basis)
    call require(size(basis) == expected_dimension, &
         'wrong adjoint-amplitude array size')
    call require_close(matrix2, expected_all_minus, 'all-minus result')
    call amplitude%evaluate(momenta, helicities, recursive_matrix2, &
         adjoint_amplitudes=recursive_basis, use_analytic_mhv=.false.)
    call require_close(recursive_matrix2, expected_all_minus, &
         'forced-BG all-minus result')
    call amplitude%evaluate(momenta, helicities, analytic_matrix2, &
         adjoint_amplitudes=analytic_basis, use_analytic_mhv=.true.)
    call require_close(analytic_matrix2, recursive_matrix2, &
         'analytic and BG all-minus results')
    call require_complex_vector_close(basis, analytic_basis, &
         'default and explicit analytic MHV vectors')
    call require_complex_vector_close(analytic_basis, recursive_basis, &
         'analytic and BG MHV vectors')

    do leg = 1, total_gluons
      if (mod(leg, 2) == 1) then
        helicities(leg) = -1
      else
        helicities(leg) = 1
      end if
    end do
    call amplitude%evaluate(momenta, helicities, matrix2)
    call require_close(matrix2, expected_alternating, 'alternating result')
    call amplitude%evaluate(momenta, helicities, recursive_matrix2, &
         use_analytic_mhv=.false.)
    call require_close(recursive_matrix2, matrix2, &
         'default and forced-BG alternating results')
  end subroutine check_multiplicity

  subroutine check_mhv_switch()
    type(ampligluon_adjoint_t) :: amplitude
    integer :: first, gluon, helicities(5), outgoing_helicities(5), second
    real(dp) :: momenta(0:3, 5)

    call fill_momenta(5, momenta)
    call amplitude%initialize(3)
    do first = 1, 4
      do second = first+1, 5
        outgoing_helicities = -1
        outgoing_helicities(first) = 1
        outgoing_helicities(second) = 1
        do gluon = 1, 5
          helicities(gluon) = outgoing_helicities(gluon)
          if (gluon <= 2) helicities(gluon) = -helicities(gluon)
        end do
        call compare_mhv_setting(amplitude, momenta, helicities, &
             'MHV switch mismatch')
        helicities = -helicities
        call compare_mhv_setting(amplitude, momenta, helicities, &
             'anti-MHV switch mismatch')
      end do
    end do

    outgoing_helicities = [-1, 1, -1, 1, -1]
    do gluon = 1, 5
      helicities(gluon) = outgoing_helicities(gluon)
      if (gluon <= 2) helicities(gluon) = -helicities(gluon)
    end do
    call compare_mhv_vector_setting(amplitude, momenta, helicities, &
         'MHV vector switch mismatch')
    helicities = -helicities
    call compare_mhv_vector_setting(amplitude, momenta, helicities, &
         'anti-MHV vector switch mismatch')
  end subroutine check_mhv_switch

  subroutine compare_mhv_vector_setting(amplitude, momenta, helicities, &
       description)
    type(ampligluon_adjoint_t), intent(inout) :: amplitude
    real(dp), intent(in) :: momenta(0:, :)
    integer, intent(in) :: helicities(:)
    character(len=*), intent(in) :: description

    complex(dp), allocatable :: analytic_basis(:), recursive_basis(:)
    real(dp) :: analytic_matrix2, recursive_matrix2

    call amplitude%evaluate(momenta, helicities, analytic_matrix2, &
         strong_coupling=0.73_dp, adjoint_amplitudes=analytic_basis, &
         use_analytic_mhv=.true.)
    call amplitude%evaluate(momenta, helicities, recursive_matrix2, &
         strong_coupling=0.73_dp, adjoint_amplitudes=recursive_basis, &
         use_analytic_mhv=.false.)
    call require_close(analytic_matrix2, recursive_matrix2, description)
    call require_complex_vector_close(analytic_basis, recursive_basis, &
         description)
  end subroutine compare_mhv_vector_setting

  subroutine compare_mhv_setting(amplitude, momenta, helicities, description)
    type(ampligluon_adjoint_t), intent(inout) :: amplitude
    real(dp), intent(in) :: momenta(0:, :)
    integer, intent(in) :: helicities(:)
    character(len=*), intent(in) :: description

    real(dp) :: analytic_matrix2, recursive_matrix2

    call amplitude%evaluate(momenta, helicities, analytic_matrix2, &
         strong_coupling=0.73_dp, average_initial_colours=.true., &
         use_analytic_mhv=.true.)
    call amplitude%evaluate(momenta, helicities, recursive_matrix2, &
         strong_coupling=0.73_dp, average_initial_colours=.true., &
         use_analytic_mhv=.false.)
    call require_close(analytic_matrix2, recursive_matrix2, description)
  end subroutine compare_mhv_setting

  subroutine check_options()
    type(ampligluon_adjoint_t) :: amplitude
    complex(dp), allocatable :: basis(:)
    real(dp) :: matrix2, momenta(0:3, 4), recursive_matrix2
    integer :: helicities(4)

    call fill_momenta(4, momenta)
    call amplitude%initialize(2)
    helicities = -1
    call amplitude%evaluate(momenta, helicities, matrix2, &
         average_initial_colours=.true.)
    call require_close(matrix2, 54.0_dp, 'initial-colour average')
    helicities = 1
    call amplitude%evaluate(momenta, helicities, matrix2)
    call require_close(matrix2, 3456.0_dp, 'analytic all-plus result')
    call amplitude%evaluate(momenta, helicities, recursive_matrix2, &
         use_analytic_mhv=.false.)
    call require_close(recursive_matrix2, matrix2, &
         'analytic and BG all-plus results')
    helicities = -1
    call amplitude%evaluate(momenta, helicities, matrix2, &
         strong_coupling=0.5_dp)
    call require_close(matrix2, 216.0_dp, 'coupling power')
    call amplitude%evaluate(momenta, helicities, matrix2, &
         strong_coupling=0.0_dp, adjoint_amplitudes=basis)
    call require(abs(matrix2) < tiny(1.0_dp), 'zero coupling is nonzero')
    call require(maxval(abs(basis)) < tiny(1.0_dp), &
         'zero-coupling adjoint amplitudes are nonzero')
    helicities = [1, 1, -1, -1]
    call amplitude%evaluate(momenta, helicities, matrix2, &
         adjoint_amplitudes=basis)
    call require(abs(matrix2) < tiny(1.0_dp), &
         'forbidden helicity is nonzero')
    call require(maxval(abs(basis)) < tiny(1.0_dp), &
         'forbidden helicity adjoint amplitudes are nonzero')
  end subroutine check_options

  subroutine fill_momenta(total_gluons, momenta)
    integer, intent(in) :: total_gluons
    real(dp), intent(out) :: momenta(0:3, total_gluons)

    real(dp) :: energy

    momenta = 0.0_dp
    momenta(:, 1) = [500.0_dp, 0.0_dp, 0.0_dp, 500.0_dp]
    momenta(:, 2) = [500.0_dp, 0.0_dp, 0.0_dp, -500.0_dp]
    select case(total_gluons)
    case(4)
      momenta(:, 3) = [500.0_dp, 300.0_dp, 400.0_dp, 0.0_dp]
      momenta(:, 4) = [500.0_dp, -300.0_dp, -400.0_dp, 0.0_dp]
    case(5)
      energy = 1000.0_dp/3.0_dp
      momenta(:, 3) = [energy, energy, 0.0_dp, 0.0_dp]
      momenta(:, 4) = [energy, -energy/2.0_dp, &
           sqrt(3.0_dp)*energy/2.0_dp, 0.0_dp]
      momenta(:, 5) = [energy, -energy/2.0_dp, &
           -sqrt(3.0_dp)*energy/2.0_dp, 0.0_dp]
    case(6)
      momenta(:, 3) = [200.0_dp, 200.0_dp, 0.0_dp, 0.0_dp]
      momenta(:, 4) = [200.0_dp, -200.0_dp, 0.0_dp, 0.0_dp]
      momenta(:, 5) = [300.0_dp, 0.0_dp, 300.0_dp, 0.0_dp]
      momenta(:, 6) = [300.0_dp, 0.0_dp, -300.0_dp, 0.0_dp]
    case(7)
      energy = 200.0_dp
      momenta(:, 3) = [energy, energy, 0.0_dp, 0.0_dp]
      momenta(:, 4) = [energy, -energy/2.0_dp, &
           sqrt(3.0_dp)*energy/2.0_dp, 0.0_dp]
      momenta(:, 5) = [energy, -energy/2.0_dp, &
           -sqrt(3.0_dp)*energy/2.0_dp, 0.0_dp]
      momenta(:, 6) = [energy, energy/sqrt(14.0_dp), &
           2.0_dp*energy/sqrt(14.0_dp), &
           3.0_dp*energy/sqrt(14.0_dp)]
      momenta(:, 7) = [energy, -energy/sqrt(14.0_dp), &
           -2.0_dp*energy/sqrt(14.0_dp), &
           -3.0_dp*energy/sqrt(14.0_dp)]
    case(8)
      energy = 1000.0_dp/6.0_dp
      momenta(:, 3) = [energy, energy, 0.0_dp, 0.0_dp]
      momenta(:, 4) = [energy, -energy, 0.0_dp, 0.0_dp]
      momenta(:, 5) = [energy, 0.0_dp, energy, 0.0_dp]
      momenta(:, 6) = [energy, 0.0_dp, -energy, 0.0_dp]
      momenta(:, 7) = [energy, energy/sqrt(14.0_dp), &
           2.0_dp*energy/sqrt(14.0_dp), &
           3.0_dp*energy/sqrt(14.0_dp)]
      momenta(:, 8) = [energy, -energy/sqrt(14.0_dp), &
           -2.0_dp*energy/sqrt(14.0_dp), &
           -3.0_dp*energy/sqrt(14.0_dp)]
    case(9)
      momenta(:, 3) = [1.33656073016398295e2_dp, &
           -4.82122784625816578e1_dp, -5.99702474874227107e1_dp, &
           1.09284452123385549e2_dp]
      momenta(:, 4) = [1.70108173084070131e2_dp, &
           1.04234173317532736e2_dp, -1.30009480500312833e2_dp, &
           3.41988690288039834e1_dp]
      momenta(:, 5) = [2.23232322123527140e1_dp, &
           -1.89496538499354372e1_dp, -1.04638307359553160e1_dp, &
           -5.45394918417156482_dp]
      momenta(:, 6) = [2.95677940418664718e2_dp, &
           -2.42855353793315089e2_dp, 1.24800823459105885e2_dp, &
           -1.13452527728897593e2_dp]
      momenta(:, 7) = [2.19471033461719117e2_dp, &
           1.98684645570178162e2_dp, 2.84833932668446863e1_dp, &
           -8.87729826659776933e1_dp]
      momenta(:, 8) = [1.13392885499946971e2_dp, &
           3.44606070664779054e1_dp, 4.11505879702098838e1_dp, &
           9.98851448030182212e1_dp]
      momenta(:, 9) = [4.53706623068480610e1_dp, &
           -2.73621398483566338e1_dp, 6.00875402753041499_dp, &
           -3.56890063761609184e1_dp]
    case(10)
      momenta(:, 3) = [2.35135902197059778e2_dp, &
           1.31018731695941909e2_dp, -1.79988114731108055e2_dp, &
           -7.56786826150768803e1_dp]
      momenta(:, 4) = [2.85393513511393053e2_dp, &
           -2.08364429070852708e2_dp, 1.14136005531997967e2_dp, &
           1.58135051438809796e2_dp]
      momenta(:, 5) = [5.20638286700562389e1_dp, &
           3.99775783286800745e1_dp, -1.68242447796891810e1_dp, &
           -2.87989630777276489e1_dp]
      momenta(:, 6) = [6.08943399536579406e1_dp, &
           4.40296801912216438e1_dp, 2.99386284100470945e1_dp, &
           2.95497280795530770e1_dp]
      momenta(:, 7) = [9.97869438452211597e1_dp, &
           -2.99586643942999142_dp, -3.58865199630482010e1_dp, &
           -9.30624340493421727e1_dp]
      momenta(:, 8) = [9.90116290745128964e1_dp, &
           -5.00823401976048146e1_dp, 1.29624658499649676e1_dp, &
           8.44217766421010367e1_dp]
      momenta(:, 9) = [9.89806591387413306e1_dp, &
           8.04544714599448980e1_dp, 1.95089773294959308e1_dp, &
           -5.42554025807367992e1_dp]
      momenta(:, 10) = [6.87331836093575674e1_dp, &
           -3.40378259679010142e1_dp, 5.61528023523394637e1_dp, &
           -2.03110738375804090e1_dp]
    case default
      call fail('unsupported test multiplicity')
    end select
  end subroutine fill_momenta

  subroutine require_close(actual, expected, description)
    real(dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: description

    if (abs(actual-expected) > max(1.0e-24_dp, 5.0e-10_dp*abs(expected))) then
      write(*, '(a,2(1x,es24.16))') trim(description), actual, expected
      call fail('matrix-element regression mismatch')
    end if
  end subroutine require_close

  subroutine require_complex_vector_close(actual, expected, description)
    complex(dp), intent(in) :: actual(:), expected(:)
    character(len=*), intent(in) :: description

    real(dp) :: error, scale

    call require(size(actual) == size(expected), &
         'complex-vector regression shape mismatch')
    error = maxval(abs(actual-expected))
    scale = maxval(abs(expected))
    if (error > max(1.0e-24_dp, 5.0e-10_dp*scale)) then
      write(*, '(a,2(1x,es24.16))') trim(description), error, scale
      call fail('complex-vector regression mismatch')
    end if
  end subroutine require_complex_vector_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) call fail(message)
  end subroutine require

end program test_ampligluon_adjoint
