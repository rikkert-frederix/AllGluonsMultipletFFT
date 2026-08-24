module gluon_kinematics
  use ampligluon_common, only: dp, fail
  implicit none
  private

  complex(dp), parameter :: imaginary_unit = (0.0_dp, 1.0_dp)
  real(dp), parameter :: inverse_sqrt_two = 1.0_dp/sqrt(2.0_dp)

  public :: aux_tensor_gluon_to_gluon
  public :: aux_tensor_gluon_terminal
  public :: external_massless_vector
  public :: gluon_aux_tensor_to_gluon
  public :: gluon_interaction_bundle
  public :: gluon_aux_tensor_terminal
  public :: massless_vector_propagator
  public :: three_gluon_terminal
  public :: terminal_vector_contract
  public :: three_gluon
  public :: two_gluon_to_aux_tensor

contains

  subroutine gluon_interaction_bundle(left, left_momentum, left_tensor, &
       has_left_tensor, right, right_momentum, right_tensor, &
       has_right_tensor, build_tensor, three_result, tensor_result, &
       tensor_gluon_result, gluon_tensor_result)
    complex(dp), intent(in) :: left(4), left_tensor(6), right(4), right_tensor(6)
    real(dp), intent(in) :: left_momentum(0:3), right_momentum(0:3)
    logical, intent(in) :: has_left_tensor, has_right_tensor, build_tensor
    complex(dp), intent(out) :: three_result(4), tensor_result(6)
    complex(dp), intent(out) :: tensor_gluon_result(4), gluon_tensor_result(4)

    call three_gluon(left, left_momentum, right, right_momentum, three_result)
    if (build_tensor) call two_gluon_to_aux_tensor(left, right, tensor_result)
    if (has_left_tensor) &
         call aux_tensor_gluon_to_gluon(left_tensor, right, tensor_gluon_result)
    if (has_right_tensor) &
         call gluon_aux_tensor_to_gluon(left, right_tensor, gluon_tensor_result)
  end subroutine gluon_interaction_bundle

  subroutine external_massless_vector(momentum, helicity, wavefunction)
    real(dp), intent(in) :: momentum(0:3)
    integer, intent(in) :: helicity
    complex(dp), intent(out) :: wavefunction(4)

    real(dp) :: effective_helicity, momentum_norm, transverse
    real(dp) :: z_over_transverse

    if (helicity /= -1 .and. helicity /= 1) &
         call fail('a gluon helicity must be -1 or +1')
    if (abs(momentum(0)) <= tiny(1.0_dp)) &
         call fail('cannot construct a polarization for zero energy')

    wavefunction = (0.0_dp, 0.0_dp)
    if (momentum(0) > 0.0_dp) then
      effective_helicity = real(helicity, dp)
      momentum_norm = momentum(0)
      transverse = sqrt(momentum(1)**2+momentum(2)**2)
      wavefunction(4) = cmplx(effective_helicity*transverse/ &
           momentum_norm*inverse_sqrt_two, 0.0_dp, dp)
      if (transverse > tiny(1.0_dp)) then
        z_over_transverse = momentum(3)/(momentum_norm*transverse)* &
             inverse_sqrt_two*effective_helicity
        wavefunction(2) = cmplx(-momentum(1)*z_over_transverse, &
             -momentum(2)/transverse*inverse_sqrt_two, dp)
        wavefunction(3) = cmplx(-momentum(2)*z_over_transverse, &
              momentum(1)/transverse*inverse_sqrt_two, dp)
      else
        wavefunction(2) = cmplx(-effective_helicity*inverse_sqrt_two, &
             0.0_dp, dp)
        wavefunction(3) = cmplx(0.0_dp, &
             sign(inverse_sqrt_two, momentum(3)), dp)
      end if
    else
      effective_helicity = real(-helicity, dp)
      momentum_norm = -momentum(0)
      transverse = sqrt(momentum(1)**2+momentum(2)**2)
      wavefunction(4) = cmplx(effective_helicity*transverse/ &
           momentum_norm*inverse_sqrt_two, 0.0_dp, dp)
      if (transverse > tiny(1.0_dp)) then
        z_over_transverse = -momentum(3)/(momentum_norm*transverse)* &
             inverse_sqrt_two*effective_helicity
        wavefunction(2) = cmplx(momentum(1)*z_over_transverse, &
             momentum(2)/transverse*inverse_sqrt_two, dp)
        wavefunction(3) = cmplx(momentum(2)*z_over_transverse, &
             -momentum(1)/transverse*inverse_sqrt_two, dp)
      else
        wavefunction(2) = cmplx(-effective_helicity*inverse_sqrt_two, &
             0.0_dp, dp)
        wavefunction(3) = cmplx(0.0_dp, &
             -sign(inverse_sqrt_two, momentum(3)), dp)
      end if
    end if
  end subroutine external_massless_vector

  subroutine three_gluon(left, left_momentum, right, right_momentum, result)
    complex(dp), intent(in) :: left(4), right(4)
    real(dp), intent(in) :: left_momentum(0:3), right_momentum(0:3)
    complex(dp), intent(out) :: result(4)

    real(dp) :: left_imag(4), left_real(4), right_imag(4), right_real(4)
    real(dp) :: left_right_imag, left_right_real, momentum_difference
    real(dp) :: polarization_imag, polarization_real
    real(dp) :: right_left_imag, right_left_real
    real(dp) :: unrotated_imag, unrotated_real
    integer :: component

    left_real = real(left, dp)
    left_imag = aimag(left)
    right_real = real(right, dp)
    right_imag = aimag(right)
    polarization_real = &
         left_real(1)*right_real(1)-left_imag(1)*right_imag(1)- &
         left_real(2)*right_real(2)+left_imag(2)*right_imag(2)- &
         left_real(3)*right_real(3)+left_imag(3)*right_imag(3)- &
         left_real(4)*right_real(4)+left_imag(4)*right_imag(4)
    polarization_imag = &
         left_real(1)*right_imag(1)+left_imag(1)*right_real(1)- &
         left_real(2)*right_imag(2)-left_imag(2)*right_real(2)- &
         left_real(3)*right_imag(3)-left_imag(3)*right_real(3)- &
         left_real(4)*right_imag(4)-left_imag(4)*right_real(4)
    left_right_real = &
         left_real(1)*right_momentum(0)-left_real(2)*right_momentum(1)- &
         left_real(3)*right_momentum(2)-left_real(4)*right_momentum(3)
    left_right_imag = &
         left_imag(1)*right_momentum(0)-left_imag(2)*right_momentum(1)- &
         left_imag(3)*right_momentum(2)-left_imag(4)*right_momentum(3)
    right_left_real = &
         right_real(1)*left_momentum(0)-right_real(2)*left_momentum(1)- &
         right_real(3)*left_momentum(2)-right_real(4)*left_momentum(3)
    right_left_imag = &
         right_imag(1)*left_momentum(0)-right_imag(2)*left_momentum(1)- &
         right_imag(3)*left_momentum(2)-right_imag(4)*left_momentum(3)
    do component = 1, 4
      momentum_difference = &
           left_momentum(component-1)-right_momentum(component-1)
      unrotated_real = polarization_real*momentum_difference+2.0_dp*( &
           left_right_real*right_real(component)- &
           left_right_imag*right_imag(component)- &
           right_left_real*left_real(component)+ &
           right_left_imag*left_imag(component))
      unrotated_imag = polarization_imag*momentum_difference+2.0_dp*( &
           left_right_real*right_imag(component)+ &
           left_right_imag*right_real(component)- &
           right_left_real*left_imag(component)- &
           right_left_imag*left_real(component))
      result(component) = cmplx(-inverse_sqrt_two*unrotated_imag, &
           inverse_sqrt_two*unrotated_real, dp)
    end do
  end subroutine three_gluon

  pure complex(dp) function three_gluon_terminal(left, left_momentum, &
       right, right_momentum, terminal) result(value)
    complex(dp), intent(in) :: left(4), right(4), terminal(4)
    real(dp), intent(in) :: left_momentum(0:3), right_momentum(0:3)

    real(dp) :: left_imag(4), left_real(4), right_imag(4), right_real(4)
    real(dp) :: terminal_imag(4), terminal_real(4)
    real(dp) :: left_right_imag, left_right_real
    real(dp) :: polarization_imag, polarization_real
    real(dp) :: right_left_imag, right_left_real
    real(dp) :: terminal_difference_imag, terminal_difference_real
    real(dp) :: terminal_left_imag, terminal_left_real
    real(dp) :: terminal_right_imag, terminal_right_real
    real(dp) :: unrotated_imag, unrotated_real
    integer :: component

    left_real = real(left, dp)
    left_imag = aimag(left)
    right_real = real(right, dp)
    right_imag = aimag(right)
    terminal_real = real(terminal, dp)
    terminal_imag = aimag(terminal)
    polarization_real = &
         left_real(1)*right_real(1)-left_imag(1)*right_imag(1)- &
         left_real(2)*right_real(2)+left_imag(2)*right_imag(2)- &
         left_real(3)*right_real(3)+left_imag(3)*right_imag(3)- &
         left_real(4)*right_real(4)+left_imag(4)*right_imag(4)
    polarization_imag = &
         left_real(1)*right_imag(1)+left_imag(1)*right_real(1)- &
         left_real(2)*right_imag(2)-left_imag(2)*right_real(2)- &
         left_real(3)*right_imag(3)-left_imag(3)*right_real(3)- &
         left_real(4)*right_imag(4)-left_imag(4)*right_real(4)
    left_right_real = &
         left_real(1)*right_momentum(0)-left_real(2)*right_momentum(1)- &
         left_real(3)*right_momentum(2)-left_real(4)*right_momentum(3)
    left_right_imag = &
         left_imag(1)*right_momentum(0)-left_imag(2)*right_momentum(1)- &
         left_imag(3)*right_momentum(2)-left_imag(4)*right_momentum(3)
    right_left_real = &
         right_real(1)*left_momentum(0)-right_real(2)*left_momentum(1)- &
         right_real(3)*left_momentum(2)-right_real(4)*left_momentum(3)
    right_left_imag = &
         right_imag(1)*left_momentum(0)-right_imag(2)*left_momentum(1)- &
         right_imag(3)*left_momentum(2)-right_imag(4)*left_momentum(3)
    terminal_difference_real = 0.0_dp
    terminal_difference_imag = 0.0_dp
    terminal_left_real = 0.0_dp
    terminal_left_imag = 0.0_dp
    terminal_right_real = 0.0_dp
    terminal_right_imag = 0.0_dp
    do component = 1, 4
      terminal_difference_real = terminal_difference_real+ &
           terminal_real(component)*(left_momentum(component-1)- &
           right_momentum(component-1))
      terminal_difference_imag = terminal_difference_imag+ &
           terminal_imag(component)*(left_momentum(component-1)- &
           right_momentum(component-1))
      terminal_left_real = terminal_left_real+ &
           terminal_real(component)*left_real(component)- &
           terminal_imag(component)*left_imag(component)
      terminal_left_imag = terminal_left_imag+ &
           terminal_real(component)*left_imag(component)+ &
           terminal_imag(component)*left_real(component)
      terminal_right_real = terminal_right_real+ &
           terminal_real(component)*right_real(component)- &
           terminal_imag(component)*right_imag(component)
      terminal_right_imag = terminal_right_imag+ &
           terminal_real(component)*right_imag(component)+ &
           terminal_imag(component)*right_real(component)
    end do
    unrotated_real = &
         polarization_real*terminal_difference_real- &
         polarization_imag*terminal_difference_imag+2.0_dp*( &
         left_right_real*terminal_right_real- &
         left_right_imag*terminal_right_imag- &
         right_left_real*terminal_left_real+ &
         right_left_imag*terminal_left_imag)
    unrotated_imag = &
         polarization_real*terminal_difference_imag+ &
         polarization_imag*terminal_difference_real+2.0_dp*( &
         left_right_real*terminal_right_imag+ &
         left_right_imag*terminal_right_real- &
         right_left_real*terminal_left_imag- &
         right_left_imag*terminal_left_real)
    value = cmplx(-inverse_sqrt_two*unrotated_imag, &
         inverse_sqrt_two*unrotated_real, dp)
  end function three_gluon_terminal

  subroutine two_gluon_to_aux_tensor(left, right, result)
    complex(dp), intent(in) :: left(4), right(4)
    complex(dp), intent(out) :: result(6)

    real(dp) :: left_imag(4), left_real(4), right_imag(4), right_real(4)

    left_real = real(left, dp)
    left_imag = aimag(left)
    right_real = real(right, dp)
    right_imag = aimag(right)
    call antisymmetric_product(1, 2, left_real, left_imag, right_real, &
         right_imag, result(1))
    call antisymmetric_product(1, 3, left_real, left_imag, right_real, &
         right_imag, result(2))
    call antisymmetric_product(1, 4, left_real, left_imag, right_real, &
         right_imag, result(3))
    call antisymmetric_product(2, 3, left_real, left_imag, right_real, &
         right_imag, result(4))
    call antisymmetric_product(2, 4, left_real, left_imag, right_real, &
         right_imag, result(5))
    call antisymmetric_product(3, 4, left_real, left_imag, right_real, &
         right_imag, result(6))
  end subroutine two_gluon_to_aux_tensor

  subroutine aux_tensor_gluon_to_gluon(left, right, result)
    complex(dp), intent(in) :: left(6), right(4)
    complex(dp), intent(out) :: result(4)

    complex(dp) :: unrotated(4)

    unrotated(1) = &
         left(1)*right(2)+left(2)*right(3)+left(3)*right(4)
    unrotated(2) = &
         left(1)*right(1)+left(4)*right(3)+left(5)*right(4)
    unrotated(3) = &
         left(2)*right(1)-left(4)*right(2)+left(6)*right(4)
    unrotated(4) = &
         left(3)*right(1)-left(5)*right(2)-left(6)*right(3)
    result = cmplx(-0.5_dp*aimag(unrotated), &
         0.5_dp*real(unrotated, dp), dp)
  end subroutine aux_tensor_gluon_to_gluon

  pure complex(dp) function aux_tensor_gluon_terminal(left, right, &
       terminal) result(value)
    complex(dp), intent(in) :: left(6), right(4), terminal(4)

    value = imaginary_unit*0.5_dp*( &
         terminal(1)*(left(1)*right(2)+left(2)*right(3)+left(3)*right(4))+ &
         terminal(2)*(left(1)*right(1)+left(4)*right(3)+left(5)*right(4))+ &
         terminal(3)*(left(2)*right(1)-left(4)*right(2)+left(6)*right(4))+ &
         terminal(4)*(left(3)*right(1)-left(5)*right(2)-left(6)*right(3)))
  end function aux_tensor_gluon_terminal

  subroutine gluon_aux_tensor_to_gluon(left, right, result)
    complex(dp), intent(in) :: left(4), right(6)
    complex(dp), intent(out) :: result(4)

    complex(dp) :: unrotated(4)

    unrotated(1) = &
         -left(2)*right(1)-left(3)*right(2)-left(4)*right(3)
    unrotated(2) = &
         -left(1)*right(1)-left(3)*right(4)-left(4)*right(5)
    unrotated(3) = &
         -left(1)*right(2)+left(2)*right(4)-left(4)*right(6)
    unrotated(4) = &
         -left(1)*right(3)+left(2)*right(5)+left(3)*right(6)
    result = cmplx(-0.5_dp*aimag(unrotated), &
         0.5_dp*real(unrotated, dp), dp)
  end subroutine gluon_aux_tensor_to_gluon

  pure complex(dp) function gluon_aux_tensor_terminal(left, right, &
       terminal) result(value)
    complex(dp), intent(in) :: left(4), right(6), terminal(4)

    value = imaginary_unit*0.5_dp*( &
         terminal(1)*(-left(2)*right(1)-left(3)*right(2)-left(4)*right(3))+ &
         terminal(2)*(-left(1)*right(1)-left(3)*right(4)-left(4)*right(5))+ &
         terminal(3)*(-left(1)*right(2)+left(2)*right(4)-left(4)*right(6))+ &
         terminal(4)*(-left(1)*right(3)+left(2)*right(5)+left(3)*right(6)))
  end function gluon_aux_tensor_terminal

  subroutine massless_vector_propagator(current, momentum)
    complex(dp), intent(inout) :: current(4)
    real(dp), intent(in) :: momentum(0:3)

    real(dp) :: virtuality

    virtuality = momentum(0)**2-momentum(1)**2- &
         momentum(2)**2-momentum(3)**2
    if (abs(virtuality) <= tiny(1.0_dp)) &
         call fail('encountered an exactly on-shell internal gluon')
    current = (-imaginary_unit/virtuality)*current
  end subroutine massless_vector_propagator

  pure complex(dp) function terminal_vector_contract(left, right) result(value)
    complex(dp), intent(in) :: left(4), right(4)

    ! AmpliCol's terminal current already carries the vector-index metric.
    value = sum(left*right)
  end function terminal_vector_contract

  pure subroutine antisymmetric_product(first, second, left_real, left_imag, &
       right_real, right_imag, result)
    integer, intent(in) :: first, second
    real(dp), intent(in) :: left_real(4), left_imag(4)
    real(dp), intent(in) :: right_real(4), right_imag(4)
    complex(dp), intent(out) :: result

    result = cmplx( &
         left_real(first)*right_real(second)- &
         left_imag(first)*right_imag(second)- &
         left_real(second)*right_real(first)+ &
         left_imag(second)*right_imag(first), &
         left_real(first)*right_imag(second)+ &
         left_imag(first)*right_real(second)- &
         left_real(second)*right_imag(first)- &
         left_imag(second)*right_real(first), dp)
  end subroutine antisymmetric_product

end module gluon_kinematics
