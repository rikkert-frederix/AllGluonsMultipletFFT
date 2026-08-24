program test_ampligluon_multiplet
  ! The nonzero reference values below were generated independently with
  ! AmpliCol imode=2 followed by its full-colour Gram-matrix contraction.
  use ampligluon_multiplet, only: ampligluon_multiplet_t
  use ampligluon_multiplet_kinds, only: dp, fail
  implicit none

  character(len=1024) :: maximum_argument, table_file
  integer :: argument_status, maximum_total_gluons

  if (command_argument_count() < 1 .or. command_argument_count() > 2) &
       call fail('usage: test_ampligluon_multiplet WIGNER_TABLE [MAX_TOTAL_GLUONS]')
  call get_command_argument(1, table_file)
  maximum_total_gluons = 8
  if (command_argument_count() == 2) then
    call get_command_argument(2, maximum_argument)
    read(maximum_argument, *, iostat=argument_status) maximum_total_gluons
    if (argument_status /= 0 .or. maximum_total_gluons < 4 .or. &
        maximum_total_gluons > 8) &
         call fail('MAX_TOTAL_GLUONS must be between four and eight')
  end if

  call check_four_gluons(trim(table_file))
  if (maximum_total_gluons >= 5) call check_five_gluons(trim(table_file))
  if (maximum_total_gluons >= 5) call check_near_collinear_mhv(trim(table_file))
  if (maximum_total_gluons >= 6) call check_six_gluons(trim(table_file))
  if (maximum_total_gluons >= 7) call check_seven_gluons(trim(table_file))
  if (maximum_total_gluons >= 8) call check_eight_gluons(trim(table_file))
  write(*, '(a)') 'AmpliGluonMultiplet regression: PASS'

contains

  subroutine check_four_gluons(filename)
    character(len=*), intent(in) :: filename

    type(ampligluon_multiplet_t) :: amplitude
    complex(dp), allocatable :: basis(:)
    real(dp) :: matrix2, momenta(0:3, 4)
    integer :: helicities(4)

    momenta(:, 1) = [500.0_dp, 0.0_dp, 0.0_dp, 500.0_dp]
    momenta(:, 2) = [500.0_dp, 0.0_dp, 0.0_dp, -500.0_dp]
    momenta(:, 3) = [500.0_dp, 300.0_dp, 400.0_dp, 0.0_dp]
    momenta(:, 4) = [500.0_dp, -300.0_dp, -400.0_dp, 0.0_dp]
    call amplitude%initialize(2, filename)
    call require(amplitude%number_of_final_gluons() == 2, &
                 'wrong final-gluon count')
    call require(amplitude%number_of_basis_amplitudes() == 8, &
                 'wrong four-gluon basis size')

    helicities = -1
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require_close(matrix2, 3456.0_dp, 'four-gluon all-minus result')
    call require_complex_close(basis(1), cmplx(0.0_dp, 24.0_dp, dp), &
                               'four-gluon basis path 1')
    call require_complex_close(basis(3), &
         cmplx(0.0_dp, 24.0_dp*sqrt(2.0_dp), dp), &
         'four-gluon basis path 3')
    call require_complex_close(basis(7), &
         cmplx(0.0_dp, -24.0_dp*sqrt(3.0_dp), dp), &
         'four-gluon basis path 7')

    call amplitude%evaluate(momenta, helicities, basis, matrix2, &
                            average_initial_colours=.true.)
    call require_close(matrix2, 54.0_dp, 'initial-colour average')

    call amplitude%evaluate(momenta, helicities, basis, matrix2, &
                            strong_coupling=0.5_dp)
    call require_close(matrix2, 216.0_dp, 'four-gluon coupling power')

    call amplitude%evaluate(momenta, helicities, basis, matrix2, &
                            strong_coupling=0.0_dp)
    call require(matrix2 <= 0.0_dp .and. maxval(abs(basis)) <= 0.0_dp, &
                 'zero coupling did not return an exact zero')

    call amplitude%evaluate(momenta, helicities, basis, matrix2, &
                            strong_coupling=-0.5_dp)
    call require_close(matrix2, 216.0_dp, 'negative coupling power')

    helicities = [1, -1, 1, -1]
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require_close(matrix2, 216.0_dp, 'four-gluon alternating helicities')

    helicities = [1, 1, -1, -1]
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require(matrix2 < 2.0e-25_dp, 'forbidden four-gluon helicity is nonzero')
    call require(maxval(abs(basis)) <= 0.0_dp, &
                 'forbidden four-gluon basis amplitude is nonzero')
  end subroutine check_four_gluons

  subroutine check_five_gluons(filename)
    character(len=*), intent(in) :: filename

    type(ampligluon_multiplet_t) :: amplitude
    complex(dp), allocatable :: basis(:)
    real(dp) :: energy, matrix2, momenta(0:3, 5), temporary(0:3)
    integer :: helicities(5)

    energy = 1000.0_dp/3.0_dp
    momenta(:, 1) = [500.0_dp, 0.0_dp, 0.0_dp, 500.0_dp]
    momenta(:, 2) = [500.0_dp, 0.0_dp, 0.0_dp, -500.0_dp]
    momenta(:, 3) = [energy, energy, 0.0_dp, 0.0_dp]
    momenta(:, 4) = [energy, -energy/2.0_dp, &
                     sqrt(3.0_dp)*energy/2.0_dp, 0.0_dp]
    momenta(:, 5) = [energy, -energy/2.0_dp, &
                    -sqrt(3.0_dp)*energy/2.0_dp, 0.0_dp]
    call amplitude%initialize(3, filename)
    call require(amplitude%number_of_basis_amplitudes() == 32, &
                 'wrong five-gluon basis size')

    helicities = -1
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require_close(matrix2, 0.839808_dp, 'five-gluon all-minus result')

    temporary = momenta(:, 3)
    momenta(:, 3) = momenta(:, 5)
    momenta(:, 5) = temporary
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require_close(matrix2, 0.839808_dp, &
                       'five-gluon closing-leg permutation')
    temporary = momenta(:, 3)
    momenta(:, 3) = momenta(:, 5)
    momenta(:, 5) = temporary

    helicities = [-1, 1, -1, 1, -1]
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require_close(matrix2, 0.010368_dp, &
                       'five-gluon alternating helicities')
  end subroutine check_five_gluons

  subroutine check_near_collinear_mhv(filename)
    character(len=*), intent(in) :: filename

    type(ampligluon_multiplet_t) :: amplitude
    complex(dp), allocatable :: basis(:)
    real(dp), parameter :: theta = 8.0e-4_dp
    real(dp) :: matrix2, momenta(0:3, 5)
    integer :: helicities(5)

    ! Legs one and three are almost collinear but nonadjacent in the planar
    ! phase calibration. The analytic inverse-soft sum is deliberately sent
    ! to the general diagonal-multiplet recursion in this ill-conditioned
    ! corner; these clean-HEAD coefficients ensure that fallback is exact.
    momenta(:, 1) = [500.0_dp, 0.0_dp, 0.0_dp, 500.0_dp]
    momenta(:, 2) = [500.0_dp, 0.0_dp, 0.0_dp, -500.0_dp]
    momenta(:, 3) = [200.0_dp, 200.0_dp*sin(theta), 0.0_dp, &
                     200.0_dp*cos(theta)]
    momenta(:, 4) = [400.0_dp, -100.0_dp*sin(theta), &
                     sqrt(600000.0_dp)/2.0_dp, -100.0_dp*cos(theta)]
    momenta(:, 5) = [400.0_dp, -100.0_dp*sin(theta), &
                    -sqrt(600000.0_dp)/2.0_dp, -100.0_dp*cos(theta)]
    helicities = [-1, -1, -1, -1, 1]
    call amplitude%initialize(3, filename)
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require_close(matrix2, 3.6000003576495569e3_dp, &
                       'near-collinear MHV fallback result')
    call require_complex_close(basis(2), &
         cmplx(-1.7708729946044895e-2_dp, -2.9263222511168078e-5_dp, dp), &
         'near-collinear MHV basis path 2')
    call require_complex_close(basis(4), &
         cmplx(-2.8621658527763900e-2_dp, 1.7320485534133415e1_dp, dp), &
         'near-collinear MHV basis path 4')
  end subroutine check_near_collinear_mhv

  subroutine check_six_gluons(filename)
    character(len=*), intent(in) :: filename

    type(ampligluon_multiplet_t) :: amplitude
    complex(dp), allocatable :: basis(:), bg_basis(:), mhv_basis(:)
    real(dp) :: bg_matrix2, matrix2, momenta(0:3, 6)
    integer :: helicities(6)

    momenta(:, 1) = [500.0_dp, 0.0_dp, 0.0_dp, 500.0_dp]
    momenta(:, 2) = [500.0_dp, 0.0_dp, 0.0_dp, -500.0_dp]
    momenta(:, 3) = [200.0_dp, 200.0_dp, 0.0_dp, 0.0_dp]
    momenta(:, 4) = [200.0_dp, -200.0_dp, 0.0_dp, 0.0_dp]
    momenta(:, 5) = [300.0_dp, 0.0_dp, 300.0_dp, 0.0_dp]
    momenta(:, 6) = [300.0_dp, 0.0_dp, -300.0_dp, 0.0_dp]
    call amplitude%initialize(4, filename)
    call require(amplitude%number_of_basis_amplitudes() == 145, &
                 'wrong six-gluon basis size')

    helicities = -1
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require_close(matrix2, 9.3333333333333397e-4_dp, &
                       'six-gluon all-minus result')
    mhv_basis = basis

    call amplitude%evaluate(momenta, helicities, bg_basis, bg_matrix2, &
                            use_mhv_optimization=.false.)
    call require_close(bg_matrix2, matrix2, &
                       'six-gluon forced-BG MHV result')
    call require_vector_close(bg_basis, mhv_basis, &
                              'six-gluon forced-BG MHV basis')

    helicities = [-1, 1, -1, 1, -1, 1]
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require_close(matrix2, 1.2049861549399821e-5_dp, &
                       'six-gluon alternating helicities')

    call amplitude%evaluate(momenta, helicities, basis, matrix2, &
                            strong_coupling=0.0_dp)
    call require(matrix2 <= 0.0_dp .and. maxval(abs(basis)) <= 0.0_dp, &
                 'generic-helicity zero coupling did not return exact zero')

    helicities = [1, 1, -1, -1, -1, -1]
    call amplitude%evaluate(momenta, helicities, bg_basis, bg_matrix2, &
                            use_mhv_optimization=.false.)
    call require(bg_matrix2 <= 0.0_dp .and. maxval(abs(bg_basis)) <= 0.0_dp, &
                 'forced-BG switch bypassed the exact-zero selection rule')

    ! Generic workspaces were initialized lazily by the earlier nonzero call.
    ! Re-enter the MHV path to catch accidental shared-workspace state.
    helicities = -1
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require_close(matrix2, 9.3333333333333397e-4_dp, &
                       'six-gluon MHV result after generic evaluation')
  end subroutine check_six_gluons

  subroutine check_seven_gluons(filename)
    character(len=*), intent(in) :: filename

    type(ampligluon_multiplet_t) :: amplitude
    complex(dp), allocatable :: basis(:)
    real(dp) :: energy, matrix2, momenta(0:3, 7)
    integer :: helicities(7)

    energy = 200.0_dp
    momenta(:, 1) = [500.0_dp, 0.0_dp, 0.0_dp, 500.0_dp]
    momenta(:, 2) = [500.0_dp, 0.0_dp, 0.0_dp, -500.0_dp]
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
    call amplitude%initialize(5, filename)
    call require(amplitude%number_of_basis_amplitudes() == 702, &
                 'wrong seven-gluon basis size')

    helicities = -1
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require_close(matrix2, 6.9826518403237469e-6_dp, &
                       'seven-gluon all-minus result')

    helicities = [-1, 1, -1, 1, -1, 1, -1]
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require_close(matrix2, 5.8326287358333135e-9_dp, &
                       'seven-gluon alternating helicities')
  end subroutine check_seven_gluons

  subroutine check_eight_gluons(filename)
    character(len=*), intent(in) :: filename

    type(ampligluon_multiplet_t) :: amplitude
    complex(dp), allocatable :: basis(:)
    real(dp) :: energy, matrix2, momenta(0:3, 8)
    integer :: helicities(8)

    energy = 1000.0_dp/6.0_dp
    momenta(:, 1) = [500.0_dp, 0.0_dp, 0.0_dp, 500.0_dp]
    momenta(:, 2) = [500.0_dp, 0.0_dp, 0.0_dp, -500.0_dp]
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
    call amplitude%initialize(6, filename)
    call require(amplitude%number_of_basis_amplitudes() == 3598, &
                 'wrong eight-gluon basis size')

    helicities = -1
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require_close(matrix2, 2.0200927180151527e-8_dp, &
                       'eight-gluon all-minus result')

    helicities = [-1, 1, -1, 1, -1, 1, -1, 1]
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require_close(matrix2, 5.8732224850953953e-10_dp, &
                       'eight-gluon alternating helicities')

    helicities = [1, 1, -1, -1, -1, -1, -1, -1]
    call amplitude%evaluate(momenta, helicities, basis, matrix2)
    call require(matrix2 < 2.0e-35_dp, &
                 'forbidden eight-gluon helicity is nonzero')
    call require(maxval(abs(basis)) <= 0.0_dp, &
                 'forbidden eight-gluon basis amplitude is nonzero')
  end subroutine check_eight_gluons

  subroutine require_close(actual, expected, description)
    real(dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: description

    real(dp), parameter :: relative_tolerance = 3.0e-11_dp

    if (abs(actual-expected) > &
        max(1.0e-24_dp, relative_tolerance*abs(expected))) then
      write(*, '(a,2(1x,es24.16))') trim(description), actual, expected
      call fail('matrix-element regression mismatch')
    end if
  end subroutine require_close

  subroutine require_complex_close(actual, expected, description)
    complex(dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: description

    if (abs(actual-expected) > 3.0e-11_dp*max(1.0_dp, abs(expected))) then
      write(*, '(a,4(1x,es24.16))') trim(description), &
           real(actual, dp), aimag(actual), real(expected, dp), aimag(expected)
      call fail('basis-amplitude convention mismatch')
    end if
  end subroutine require_complex_close

  subroutine require_vector_close(actual, expected, description)
    complex(dp), intent(in) :: actual(:), expected(:)
    character(len=*), intent(in) :: description

    if (size(actual) /= size(expected)) call fail('basis-vector size mismatch')
    if (maxval(abs(actual-expected)) > &
        3.0e-11_dp*max(1.0_dp, maxval(abs(expected)))) then
      write(*, '(a,1x,es24.16)') trim(description), &
           maxval(abs(actual-expected))
      call fail('basis-vector regression mismatch')
    end if
  end subroutine require_vector_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) call fail(message)
  end subroutine require

end program test_ampligluon_multiplet
