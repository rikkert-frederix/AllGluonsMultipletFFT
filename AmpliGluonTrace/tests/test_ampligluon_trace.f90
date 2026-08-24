program test_ampligluon_trace
  use ampligluon_common, only: dp, fail
  use ampligluon_trace, only: ampligluon_trace_t
  use trace_colour_kernel, only: build_u3_trace_kernel
  implicit none

  call check_four_gluons()
  call check_five_gluons()
  call check_five_gluon_mhv_switch()
  call check_six_gluons()
  call check_colour_contraction_switch()
  call check_seven_gluons()
  call check_eight_gluons()
  call check_nine_gluons()
  call check_eleven_gluon_zeros()
  call check_reinitialization()
  write(*, '(a)') 'AmpliGluonTrace regression: PASS'

contains

  subroutine check_four_gluons()
    type(ampligluon_trace_t) :: amplitude
    complex(dp), allocatable :: ordered(:)
    real(dp) :: matrix2, momenta(0:3, 4)
    integer :: helicities(4)

    call fill_momenta(4, momenta)
    call amplitude%initialize(2)
    call require(amplitude%number_of_final_gluons() == 2, &
         'wrong four-gluon final-state count')
    call require(amplitude%number_of_colour_orders() == 6, &
         'wrong four-gluon colour-order count')

    helicities = -1
    call amplitude%evaluate(momenta, helicities, matrix2, &
         ordered_amplitudes=ordered)
    call require(size(ordered) == 6, 'wrong ordered-amplitude array size')
    call require_close(matrix2, 3456.0_dp, 'four-gluon all-minus result')
    call check_ordered_identities(4, ordered)
    call check_u3_colour_sum(4, ordered, matrix2)

    call amplitude%evaluate(momenta, helicities, matrix2, &
         average_initial_colours=.true.)
    call require_close(matrix2, 54.0_dp, 'initial-colour average')

    call amplitude%evaluate(momenta, helicities, matrix2, &
         strong_coupling=0.5_dp)
    call require_close(matrix2, 216.0_dp, 'four-gluon coupling power')

    helicities = [1, -1, 1, -1]
    call amplitude%evaluate(momenta, helicities, matrix2)
    call require_close(matrix2, 216.0_dp, &
         'four-gluon alternating-helicity result')

    helicities = [1, 1, -1, -1]
    call amplitude%evaluate(momenta, helicities, matrix2, &
         ordered_amplitudes=ordered)
    call require(abs(matrix2) <= 0.0_dp, &
         'forbidden four-gluon helicity is not exactly zero')
    call require(all(abs(ordered) <= 0.0_dp), &
         'forbidden four-gluon ordered amplitudes are not exactly zero')

    helicities = [-1, 1, -1, -1]
    call amplitude%evaluate(momenta, helicities, matrix2)
    call require(abs(matrix2) <= 0.0_dp, &
         'single-positive all-outgoing helicity is not exactly zero')

    helicities = [1, -1, 1, 1]
    call amplitude%evaluate(momenta, helicities, matrix2)
    call require(abs(matrix2) <= 0.0_dp, &
         'single-negative all-outgoing helicity is not exactly zero')

    helicities = [-1, -1, 1, 1]
    call amplitude%evaluate(momenta, helicities, matrix2)
    call require(abs(matrix2) <= 0.0_dp, &
         'all-positive all-outgoing helicity is not exactly zero')
  end subroutine check_four_gluons

  subroutine check_five_gluons()
    type(ampligluon_trace_t) :: amplitude
    complex(dp), allocatable :: bg_ordered(:), ordered(:), ordered_reference(:)
    real(dp) :: matrix2, matrix2_reference, momenta(0:3, 5), scale
    integer :: helicities(5)

    call fill_momenta(5, momenta)
    call amplitude%initialize(3)
    call require(amplitude%number_of_colour_orders() == 24, &
         'wrong five-gluon colour-order count')

    helicities = -1
    call amplitude%evaluate(momenta, helicities, matrix2, &
         ordered_amplitudes=ordered)
    call require_close(matrix2, 0.839808_dp, 'five-gluon all-minus result')
    call check_ordered_identities(5, ordered)
    call check_u3_colour_sum(5, ordered, matrix2)
    matrix2_reference = matrix2
    ordered_reference = ordered
    call amplitude%evaluate(momenta, helicities, matrix2, &
         strong_coupling=-0.5_dp, ordered_amplitudes=ordered)
    call require_close(matrix2, matrix2_reference*0.5_dp**6, &
         'five-gluon negative-coupling result')
    scale = max(maxval(abs(ordered)), maxval(abs(ordered_reference)), &
         tiny(1.0_dp))
    call require(maxval(abs(ordered-(-0.5_dp)**3*ordered_reference)) <= &
         5.0e-11_dp*scale, 'five-gluon ordered coupling scaling failed')
    call amplitude%evaluate(momenta, helicities, matrix2, &
         strong_coupling=-0.5_dp, ordered_amplitudes=bg_ordered, &
         use_mhv_optimization=.false.)
    call require_close(matrix2, matrix2_reference*0.5_dp**6, &
         'five-gluon forced-BG negative-coupling result')
    scale = max(maxval(abs(bg_ordered)), maxval(abs(ordered_reference)), &
         tiny(1.0_dp))
    call require(maxval(abs(bg_ordered- &
         (-0.5_dp)**3*ordered_reference)) <= 5.0e-11_dp*scale, &
         'five-gluon forced-BG ordered coupling scaling failed')

    helicities = [-1, 1, -1, 1, -1]
    call amplitude%evaluate(momenta, helicities, matrix2)
    call require_close(matrix2, 0.010368_dp, &
         'five-gluon alternating-helicity result')
  end subroutine check_five_gluons

  subroutine check_five_gluon_mhv_switch()
    type(ampligluon_trace_t) :: amplitude
    real(dp) :: momenta(0:3, 5)
    integer :: first, gluon, helicities(5), outgoing_helicities(5), second

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
             'five-gluon MHV switch mismatch')
        call compare_mhv_setting(amplitude, momenta, -helicities, &
             'five-gluon anti-MHV switch mismatch')
      end do
    end do
  end subroutine check_five_gluon_mhv_switch

  subroutine compare_mhv_setting(amplitude, momenta, helicities, description)
    type(ampligluon_trace_t), intent(inout) :: amplitude
    real(dp), intent(in) :: momenta(0:, :)
    integer, intent(in) :: helicities(:)
    character(len=*), intent(in) :: description

    complex(dp), allocatable :: analytic_ordered(:), bg_ordered(:)
    real(dp) :: analytic_matrix2, bg_matrix2, scale

    call amplitude%evaluate(momenta, helicities, analytic_matrix2, &
         strong_coupling=0.73_dp, average_initial_colours=.true., &
         ordered_amplitudes=analytic_ordered, use_mhv_optimization=.true.)
    call amplitude%evaluate(momenta, helicities, bg_matrix2, &
         strong_coupling=0.73_dp, average_initial_colours=.true., &
         ordered_amplitudes=bg_ordered, use_mhv_optimization=.false.)
    call require_close(analytic_matrix2, bg_matrix2, description)
    scale = max(maxval(abs(analytic_ordered)), maxval(abs(bg_ordered)), &
         tiny(1.0_dp))
    call require(maxval(abs(analytic_ordered-bg_ordered)) <= &
         5.0e-11_dp*scale, trim(description)//' ordered amplitudes')
  end subroutine compare_mhv_setting

  subroutine check_six_gluons()
    type(ampligluon_trace_t) :: amplitude
    complex(dp), allocatable :: analytic_ordered(:), bg_ordered(:), ordered(:)
    real(dp) :: analytic_matrix2, bg_matrix2, matrix2
    real(dp) :: momenta(0:3, 6), scale
    integer :: helicities(6)

    call fill_momenta(6, momenta)
    call amplitude%initialize(4)
    call require(amplitude%number_of_colour_orders() == 120, &
         'wrong six-gluon colour-order count')

    helicities = -1
    call amplitude%evaluate(momenta, helicities, matrix2, &
         ordered_amplitudes=ordered)
    call require_close(matrix2, 9.3333333333333397e-4_dp, &
         'six-gluon all-minus result')
    call check_ordered_identities(6, ordered)
    call check_u3_colour_sum(6, ordered, matrix2)

    call amplitude%evaluate(momenta, helicities, bg_matrix2, &
         ordered_amplitudes=bg_ordered, use_mhv_optimization=.false.)
    call require_close(bg_matrix2, matrix2, &
         'six-gluon forced-BG MHV result')
    scale = max(maxval(abs(bg_ordered)), maxval(abs(ordered)), tiny(1.0_dp))
    call require(maxval(abs(bg_ordered-ordered)) <= 5.0e-11_dp*scale, &
         'six-gluon forced-BG ordered MHV amplitudes disagree')

    helicities = 1
    call amplitude%evaluate(momenta, helicities, analytic_matrix2, &
         ordered_amplitudes=analytic_ordered)
    call amplitude%evaluate(momenta, helicities, bg_matrix2, &
         ordered_amplitudes=bg_ordered, use_mhv_optimization=.false.)
    call require_close(bg_matrix2, analytic_matrix2, &
         'six-gluon forced-BG anti-MHV result')
    scale = max(maxval(abs(bg_ordered)), maxval(abs(analytic_ordered)), &
         tiny(1.0_dp))
    call require(maxval(abs(bg_ordered-analytic_ordered)) <= &
         5.0e-11_dp*scale, &
         'six-gluon forced-BG ordered anti-MHV amplitudes disagree')

    call amplitude%evaluate(momenta, helicities, matrix2, &
         ordered_amplitudes=ordered)
    call require_close(matrix2, analytic_matrix2, &
         'six-gluon anti-MHV path after forced BG')
    scale = max(maxval(abs(ordered)), maxval(abs(analytic_ordered)), &
         tiny(1.0_dp))
    call require(maxval(abs(ordered-analytic_ordered)) <= &
         5.0e-11_dp*scale, &
         'six-gluon ordered anti-MHV path after forced BG disagrees')

    helicities = [-1, 1, -1, 1, -1, 1]
    call amplitude%evaluate(momenta, helicities, matrix2)
    call require_close(matrix2, 1.2049861549399821e-5_dp, &
         'six-gluon alternating-helicity result')
    call amplitude%evaluate(momenta, helicities, matrix2, &
         strong_coupling=-0.5_dp)
    call require_close(matrix2, 1.2049861549399821e-5_dp*0.5_dp**8, &
         'six-gluon non-MHV coupling result')
  end subroutine check_six_gluons

  subroutine check_colour_contraction_switch()
    type(ampligluon_trace_t) :: direct_amplitude, fourier_amplitude
    real(dp) :: direct_matrix2, fourier_matrix2, momenta(0:3, 6)
    integer :: helicities(6)

    call fill_momenta(6, momenta)
    call fourier_amplitude%initialize(4)
    call direct_amplitude%initialize(4, use_colour_fft=.false.)

    helicities = -1
    call fourier_amplitude%evaluate(momenta, helicities, fourier_matrix2)
    call direct_amplitude%evaluate(momenta, helicities, direct_matrix2)
    call require_close(direct_matrix2, fourier_matrix2, &
         'normal and FFT colour contractions disagree for MHV amplitudes')

    helicities = [-1, 1, -1, 1, -1, 1]
    call fourier_amplitude%evaluate(momenta, helicities, fourier_matrix2)
    call direct_amplitude%evaluate(momenta, helicities, direct_matrix2)
    call require_close(direct_matrix2, fourier_matrix2, &
         'normal and FFT colour contractions disagree for NMHV amplitudes')
  end subroutine check_colour_contraction_switch

  subroutine check_seven_gluons()
    type(ampligluon_trace_t) :: amplitude
    complex(dp), allocatable :: ordered(:)
    real(dp) :: matrix2, momenta(0:3, 7)
    integer :: helicities(7)

    call fill_momenta(7, momenta)
    call amplitude%initialize(5)
    call require(amplitude%number_of_colour_orders() == 720, &
         'wrong seven-gluon colour-order count')

    helicities = -1
    call amplitude%evaluate(momenta, helicities, matrix2, &
         ordered_amplitudes=ordered)
    call require_close(matrix2, 6.9826518403237469e-6_dp, &
         'seven-gluon all-minus result')
    call check_ordered_identities(7, ordered)
    call check_u3_colour_sum(7, ordered, matrix2)

    helicities = [-1, 1, -1, 1, -1, 1, -1]
    call amplitude%evaluate(momenta, helicities, matrix2)
    call require_close(matrix2, 5.8326287358333135e-9_dp, &
         'seven-gluon alternating-helicity result')
  end subroutine check_seven_gluons

  subroutine check_eight_gluons()
    type(ampligluon_trace_t) :: amplitude
    complex(dp), allocatable :: ordered(:)
    real(dp) :: matrix2, momenta(0:3, 8)
    integer :: helicities(8)

    call fill_momenta(8, momenta)
    call amplitude%initialize(6)
    call require(amplitude%number_of_colour_orders() == 5040, &
         'wrong eight-gluon colour-order count')

    helicities = -1
    call amplitude%evaluate(momenta, helicities, matrix2, &
         ordered_amplitudes=ordered)
    call require_close(matrix2, 2.0200927180151527e-8_dp, &
         'eight-gluon all-minus result')
    call check_ordered_identities(8, ordered)

    helicities = [-1, 1, -1, 1, -1, 1, -1, 1]
    call amplitude%evaluate(momenta, helicities, matrix2)
    call require_close(matrix2, 5.8732224850953953e-10_dp, &
         'eight-gluon alternating-helicity result')

    helicities = [1, 1, -1, -1, -1, -1, -1, -1]
    call amplitude%evaluate(momenta, helicities, matrix2, &
         ordered_amplitudes=ordered)
    call require(abs(matrix2) <= 0.0_dp, &
         'forbidden eight-gluon helicity is not exactly zero')
    call require(all(abs(ordered) <= 0.0_dp), &
         'forbidden eight-gluon ordered amplitudes are not exactly zero')
  end subroutine check_eight_gluons

  subroutine check_nine_gluons()
    type(ampligluon_trace_t) :: amplitude
    complex(dp), allocatable :: ordered(:)
    real(dp) :: matrix2, momenta(0:3, 9)
    integer :: helicities(9)

    call fill_momenta(9, momenta)
    call amplitude%initialize(7)
    call require(amplitude%number_of_colour_orders() == 40320, &
         'wrong nine-gluon colour-order count')

    helicities = -1
    call amplitude%evaluate(momenta, helicities, matrix2, &
         ordered_amplitudes=ordered)
    call require_close(matrix2, 2.9947601653227329e-8_dp, &
         'nine-gluon all-minus result')
    call check_ordered_identities(9, ordered)

    helicities = [-1, 1, -1, 1, -1, 1, -1, 1, -1]
    call amplitude%evaluate(momenta, helicities, matrix2)
    call require_close(matrix2, 1.2609237070482274e-11_dp, &
         'nine-gluon alternating-helicity result')

    helicities = [1, 1, -1, -1, -1, -1, -1, -1, -1]
    call amplitude%evaluate(momenta, helicities, matrix2, &
         ordered_amplitudes=ordered)
    call require(abs(matrix2) <= 0.0_dp, &
         'forbidden nine-gluon helicity is not exactly zero')
    call require(all(abs(ordered) <= 0.0_dp), &
         'forbidden nine-gluon ordered amplitudes are not exactly zero')
  end subroutine check_nine_gluons

  subroutine check_eleven_gluon_zeros()
    type(ampligluon_trace_t) :: amplitude
    real(dp) :: angle, energy, matrix2, momenta(0:3, 11)
    integer :: final_gluon, helicities(11)

    ! A regular massless nine-particle final state lets the high-multiplicity
    ! exact-zero path be tested without constructing any factorial colour or
    ! current data.  The spatial momenta cancel around the regular polygon.
    momenta = 0.0_dp
    momenta(:, 1) = [500.0_dp, 0.0_dp, 0.0_dp, 500.0_dp]
    momenta(:, 2) = [500.0_dp, 0.0_dp, 0.0_dp, -500.0_dp]
    energy = 1000.0_dp/9.0_dp
    do final_gluon = 1, 9
      angle = 2.0_dp*acos(-1.0_dp)*real(final_gluon-1, dp)/9.0_dp
      momenta(:, final_gluon+2) = [energy, energy*cos(angle), &
           energy*sin(angle), 0.0_dp]
    end do

    call amplitude%initialize(9, use_colour_fft=.false.)
    call require(amplitude%number_of_colour_orders() == 3628800, &
         'wrong eleven-gluon colour-order count')

    helicities = -1
    helicities(1:2) = 1
    call amplitude%evaluate(momenta, helicities, matrix2, &
         use_mhv_optimization=.false.)
    call require(abs(matrix2) <= 0.0_dp, &
         'all-negative eleven-gluon sector is not exactly zero')

    helicities = 1
    helicities(1:2) = -1
    call amplitude%evaluate(momenta, helicities, matrix2)
    call require(abs(matrix2) <= 0.0_dp, &
         'all-positive eleven-gluon sector is not exactly zero')

    helicities = [(merge(-1, 1, modulo(final_gluon, 2) == 1), &
         final_gluon=1, 11)]
    call amplitude%evaluate(momenta, helicities, matrix2, &
         strong_coupling=0.0_dp)
    call require(abs(matrix2) <= 0.0_dp, &
         'zero-coupling eleven-gluon result is not exactly zero')
  end subroutine check_eleven_gluon_zeros

  subroutine check_reinitialization()
    type(ampligluon_trace_t) :: amplitude
    complex(dp), allocatable :: ordered(:)
    real(dp) :: matrix2, momenta9(0:3, 9), momenta4(0:3, 4)
    integer :: helicities9(9), helicities4(4)

    call fill_momenta(9, momenta9)
    call amplitude%initialize(7)
    helicities9 = [-1, 1, -1, 1, -1, 1, -1, 1, -1]
    call amplitude%evaluate(momenta9, helicities9, matrix2, &
         ordered_amplitudes=ordered)
    call require_close(matrix2, 1.2609237070482274e-11_dp, &
         'nine-gluon result before reinitialization')
    call require(size(ordered) == 40320, &
         'wrong ordered size before reinitialization')

    call fill_momenta(4, momenta4)
    call amplitude%initialize(2, use_colour_fft=.false.)
    helicities4 = -1
    call amplitude%evaluate(momenta4, helicities4, matrix2, &
         ordered_amplitudes=ordered)
    call require_close(matrix2, 3456.0_dp, &
         'four-gluon result after reinitialization')
    call require(size(ordered) == 6, &
         'wrong ordered size after reinitialization')
  end subroutine check_reinitialization

  subroutine check_ordered_identities(total_gluons, ordered)
    integer, intent(in) :: total_gluons
    complex(dp), intent(in) :: ordered(:)

    integer, allocatable :: base_order(:), factorials(:), permutation(:)
    integer, allocatable :: reflected(:), with_photon(:)
    complex(dp) :: decoupling_sum, reflection_residual
    real(dp) :: decoupling_scale, reflection_scale, reflection_sign
    integer :: insertion, order, position

    allocate(factorials(0:total_gluons-1))
    factorials(0) = 1
    do position = 1, total_gluons-1
      factorials(position) = position*factorials(position-1)
    end do
    call require(size(ordered) == factorials(total_gluons-1), &
                 'wrong ordered-amplitude vector size')

    allocate(permutation(total_gluons-1), reflected(total_gluons-1))
    permutation = [(position, position=1, total_gluons-1)]
    reflection_sign = 1.0_dp
    if (modulo(total_gluons, 2) /= 0) reflection_sign = -1.0_dp
    do order = 1, size(ordered)
      reflected = permutation(size(permutation):1:-1)
      reflection_residual = ordered(order)-reflection_sign* &
           ordered(permutation_rank(reflected, factorials))
      reflection_scale = max(abs(ordered(order)), &
           abs(ordered(permutation_rank(reflected, factorials))))
      call require(abs(reflection_residual) <= &
           5.0e-11_dp*max(reflection_scale, tiny(1.0_dp)), &
           'ordered-amplitude reflection identity failed')
      if (order < size(ordered)) call next_permutation(permutation)
    end do

    ! Turn leg 1 into a U(1) gluon and insert it in every position of one
    ! fixed ordering of the other legs.  Photon decoupling makes the sum zero.
    allocate(base_order(total_gluons-2), with_photon(total_gluons-1))
    base_order = [(position, position=2, total_gluons-1)]
    decoupling_sum = (0.0_dp, 0.0_dp)
    decoupling_scale = 0.0_dp
    do insertion = 1, total_gluons-1
      with_photon(insertion) = 1
      if (insertion > 1) &
           with_photon(1:insertion-1) = base_order(1:insertion-1)
      if (insertion < total_gluons-1) &
           with_photon(insertion+1:) = base_order(insertion:)
      order = permutation_rank(with_photon, factorials)
      decoupling_sum = decoupling_sum+ordered(order)
      decoupling_scale = decoupling_scale+abs(ordered(order))
    end do
    call require(abs(decoupling_sum) <= &
         5.0e-11_dp*max(decoupling_scale, tiny(1.0_dp)), &
         'ordered-amplitude photon-decoupling identity failed')
  end subroutine check_ordered_identities

  subroutine check_u3_colour_sum(total_gluons, ordered, su3_result)
    integer, intent(in) :: total_gluons
    complex(dp), intent(in) :: ordered(:)
    real(dp), intent(in) :: su3_result

    integer, allocatable :: factorials(:), inverse(:), permutations(:, :)
    integer, allocatable :: relative(:)
    real(dp), allocatable :: kernel(:)
    real(dp) :: u3_result
    integer :: column, order, position, rank, row

    call build_u3_trace_kernel(total_gluons, kernel)
    call require(size(kernel) == size(ordered), &
                 'U(3) kernel and amplitude dimensions disagree')
    allocate(factorials(0:total_gluons-1))
    factorials(0) = 1
    do position = 1, total_gluons-1
      factorials(position) = position*factorials(position-1)
    end do
    allocate(permutations(total_gluons-1, size(ordered)))
    permutations(:, 1) = [(position, position=1, total_gluons-1)]
    do order = 2, size(ordered)
      permutations(:, order) = permutations(:, order-1)
      call next_permutation(permutations(:, order))
    end do
    allocate(inverse(total_gluons-1), relative(total_gluons-1))

    u3_result = 0.0_dp
    do row = 1, size(ordered)
      do position = 1, total_gluons-1
        inverse(permutations(position, row)) = position
      end do
      do column = 1, size(ordered)
        do position = 1, total_gluons-1
          relative(position) = inverse(permutations(position, column))
        end do
        rank = permutation_rank(relative, factorials)
        u3_result = u3_result+kernel(rank)* &
             real(ordered(column)*conjg(ordered(row)), dp)
      end do
    end do
    call require_close(u3_result, su3_result, &
                       'SU(3)--U(3) physical colour sum')
  end subroutine check_u3_colour_sum

  integer function permutation_rank(permutation, factorials) result(rank)
    integer, intent(in) :: permutation(:)
    integer, intent(in) :: factorials(0:)

    integer :: left, right, smaller

    rank = 1
    do left = 1, size(permutation)-1
      smaller = 0
      do right = left+1, size(permutation)
        if (permutation(right) < permutation(left)) smaller = smaller+1
      end do
      rank = rank+smaller*factorials(size(permutation)-left)
    end do
  end function permutation_rank

  subroutine next_permutation(permutation)
    integer, intent(inout) :: permutation(:)

    integer :: left, right, temporary

    left = size(permutation)-1
    do while (left >= 1)
      if (permutation(left) < permutation(left+1)) exit
      left = left-1
    end do
    if (left < 1) call fail('test permutation wrapped unexpectedly')
    right = size(permutation)
    do while (permutation(right) < permutation(left))
      right = right-1
    end do
    temporary = permutation(left)
    permutation(left) = permutation(right)
    permutation(right) = temporary
    permutation(left+1:) = permutation(size(permutation):left+1:-1)
  end subroutine next_permutation

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
    case default
      call fail('unsupported test multiplicity')
    end select
  end subroutine fill_momenta

  subroutine require_close(actual, expected, description)
    real(dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: description

    real(dp), parameter :: relative_tolerance = 5.0e-10_dp

    if (abs(actual-expected) > &
         max(1.0e-24_dp, relative_tolerance*abs(expected))) then
      write(*, '(a,2(1x,es24.16))') trim(description), actual, expected
      call fail('matrix-element regression mismatch')
    end if
  end subroutine require_close

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) call fail(message)
  end subroutine require

end program test_ampligluon_trace
