module qcd_kinematics
  use ampligluon_multiplet_kinds, only: dp, fail, imaginary_unit
  implicit none
  private

  real(dp), parameter :: inverse_sqrt_two = 1.0_dp/sqrt(2.0_dp)

  public :: aux_tensor_gluon_to_gluon
  public :: external_massless_vector
  public :: gluon_aux_tensor_to_gluon
  public :: lorentz_contract
  public :: massless_vector_propagator
  public :: massless_vector_propagator_copy
  public :: three_gluon
  public :: two_gluon_to_aux_tensor

contains

  subroutine external_massless_vector(momentum, helicity, wavefunction)
    real(dp), intent(in) :: momentum(0:3)
    integer, intent(in) :: helicity
    complex(dp), intent(out) :: wavefunction(4)

    real(dp) :: effective_helicity, momentum_norm, transverse, z_over_transverse

    if (helicity /= -1 .and. helicity /= 1) &
         call fail('a gluon helicity must be -1 or +1')
    if (abs(momentum(0)) <= tiny(1.0_dp)) &
         call fail('cannot construct a polarization for zero energy')

    wavefunction = (0.0_dp, 0.0_dp)
    if (momentum(0) > 0.0_dp) then
      effective_helicity = real(helicity, dp)
      momentum_norm = momentum(0)
      transverse = sqrt(momentum(1)**2 + momentum(2)**2)
      wavefunction(4) = cmplx(effective_helicity*transverse/momentum_norm* &
                              inverse_sqrt_two, 0.0_dp, dp)
      if (transverse > tiny(1.0_dp)) then
        z_over_transverse = momentum(3)/(momentum_norm*transverse)* &
                            inverse_sqrt_two*effective_helicity
        wavefunction(2) = cmplx(-momentum(1)*z_over_transverse, &
                                -momentum(2)/transverse*inverse_sqrt_two, dp)
        wavefunction(3) = cmplx(-momentum(2)*z_over_transverse, &
                                 momentum(1)/transverse*inverse_sqrt_two, dp)
      else
        wavefunction(2) = cmplx(-effective_helicity*inverse_sqrt_two, 0.0_dp, dp)
        wavefunction(3) = cmplx(0.0_dp, &
                                sign(inverse_sqrt_two, momentum(3)), dp)
      end if
    else
      effective_helicity = real(-helicity, dp)
      momentum_norm = -momentum(0)
      transverse = sqrt(momentum(1)**2 + momentum(2)**2)
      wavefunction(4) = cmplx(effective_helicity*transverse/momentum_norm* &
                              inverse_sqrt_two, 0.0_dp, dp)
      if (transverse > tiny(1.0_dp)) then
        z_over_transverse = -momentum(3)/(momentum_norm*transverse)* &
                            inverse_sqrt_two*effective_helicity
        wavefunction(2) = cmplx(momentum(1)*z_over_transverse, &
                                momentum(2)/transverse*inverse_sqrt_two, dp)
        wavefunction(3) = cmplx(momentum(2)*z_over_transverse, &
                               -momentum(1)/transverse*inverse_sqrt_two, dp)
      else
        wavefunction(2) = cmplx(-effective_helicity*inverse_sqrt_two, 0.0_dp, dp)
        wavefunction(3) = cmplx(0.0_dp, &
                               -sign(inverse_sqrt_two, momentum(3)), dp)
      end if
    end if
  end subroutine external_massless_vector

  subroutine three_gluon(left, left_momentum, right, right_momentum, result)
    complex(dp), intent(in) :: left(4), right(4)
    real(dp), intent(in) :: left_momentum(0:3), right_momentum(0:3)
    complex(dp), intent(out) :: result(4)

    complex(dp) :: polarization_product, left_with_right_momentum
    complex(dp) :: right_with_left_momentum

    polarization_product = lorentz_contract(left, right)
    left_with_right_momentum = left(1)*right_momentum(0) - &
         left(2)*right_momentum(1) - left(3)*right_momentum(2) - &
         left(4)*right_momentum(3)
    right_with_left_momentum = right(1)*left_momentum(0) - &
         right(2)*left_momentum(1) - right(3)*left_momentum(2) - &
         right(4)*left_momentum(3)
    result = imaginary_unit*inverse_sqrt_two*( &
         polarization_product*(left_momentum-right_momentum) + &
         2.0_dp*(left_with_right_momentum*right - &
                 right_with_left_momentum*left))
  end subroutine three_gluon

  subroutine two_gluon_to_aux_tensor(left, right, result)
    complex(dp), intent(in) :: left(4), right(4)
    complex(dp), intent(out) :: result(6)

    result(1) = left(1)*right(2)-left(2)*right(1)
    result(2) = left(1)*right(3)-left(3)*right(1)
    result(3) = left(1)*right(4)-left(4)*right(1)
    result(4) = left(2)*right(3)-left(3)*right(2)
    result(5) = left(2)*right(4)-left(4)*right(2)
    result(6) = left(3)*right(4)-left(4)*right(3)
  end subroutine two_gluon_to_aux_tensor

  subroutine aux_tensor_gluon_to_gluon(left, right, result)
    complex(dp), intent(in) :: left(6), right(4)
    complex(dp), intent(out) :: result(4)

    result(1) = imaginary_unit*0.5_dp*( &
         left(1)*right(2)+left(2)*right(3)+left(3)*right(4))
    result(2) = imaginary_unit*0.5_dp*( &
         left(1)*right(1)+left(4)*right(3)+left(5)*right(4))
    result(3) = imaginary_unit*0.5_dp*( &
         left(2)*right(1)-left(4)*right(2)+left(6)*right(4))
    result(4) = imaginary_unit*0.5_dp*( &
         left(3)*right(1)-left(5)*right(2)-left(6)*right(3))
  end subroutine aux_tensor_gluon_to_gluon

  subroutine gluon_aux_tensor_to_gluon(left, right, result)
    complex(dp), intent(in) :: left(4), right(6)
    complex(dp), intent(out) :: result(4)

    result(1) = imaginary_unit*0.5_dp*( &
         -left(2)*right(1)-left(3)*right(2)-left(4)*right(3))
    result(2) = imaginary_unit*0.5_dp*( &
         -left(1)*right(1)-left(3)*right(4)-left(4)*right(5))
    result(3) = imaginary_unit*0.5_dp*( &
         -left(1)*right(2)+left(2)*right(4)-left(4)*right(6))
    result(4) = imaginary_unit*0.5_dp*( &
         -left(1)*right(3)+left(2)*right(5)+left(3)*right(6))
  end subroutine gluon_aux_tensor_to_gluon

  subroutine massless_vector_propagator(current, momentum)
    complex(dp), intent(inout) :: current(:, :)
    real(dp), intent(in) :: momentum(0:3)

    complex(dp) :: propagator
    real(dp) :: virtuality

    if (size(current, 1) /= 4) &
         call fail('massless-vector current has the wrong dimension')
    virtuality = momentum(0)**2 - momentum(1)**2 - &
                 momentum(2)**2 - momentum(3)**2
    if (abs(virtuality) <= tiny(virtuality)) &
         call fail('encountered an exactly on-shell internal gluon')
    propagator = -imaginary_unit/virtuality
    current = propagator*current
  end subroutine massless_vector_propagator

  subroutine massless_vector_propagator_copy(source, current, momentum)
    complex(dp), intent(in) :: source(:, :)
    complex(dp), intent(out) :: current(:, :)
    real(dp), intent(in) :: momentum(0:3)

    complex(dp) :: propagator
    real(dp) :: virtuality

    if (size(source, 1) /= 4 .or. any(shape(current) /= shape(source))) &
         call fail('massless-vector current has the wrong dimension')
    virtuality = momentum(0)**2 - momentum(1)**2 - &
                 momentum(2)**2 - momentum(3)**2
    if (abs(virtuality) <= tiny(virtuality)) &
         call fail('encountered an exactly on-shell internal gluon')
    propagator = -imaginary_unit/virtuality
    current = propagator*source
  end subroutine massless_vector_propagator_copy

  pure complex(dp) function lorentz_contract(left, right) result(value)
    complex(dp), intent(in) :: left(4), right(4)

    value = left(1)*right(1)-left(2)*right(2)- &
            left(3)*right(3)-left(4)*right(4)
  end function lorentz_contract

end module qcd_kinematics
