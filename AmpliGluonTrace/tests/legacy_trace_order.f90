module legacy_trace_order
  use ampligluon_common, only: dp
  use gluon_kinematics, only: aux_tensor_gluon_to_gluon, &
       gluon_aux_tensor_to_gluon, massless_vector_propagator, &
       terminal_vector_contract, three_gluon, two_gluon_to_aux_tensor
  implicit none
  private

  public :: evaluate_legacy_trace_orders

contains

  subroutine evaluate_legacy_trace_orders(outgoing_momenta, &
       external_wavefunctions, coupling, amplitudes)
    real(dp), intent(in) :: outgoing_momenta(0:, :)
    complex(dp), intent(in) :: external_wavefunctions(:, :)
    real(dp), intent(in) :: coupling
    complex(dp), intent(out) :: amplitudes(:)

    integer, allocatable :: permutation(:)
    integer :: degree, order

    degree = size(outgoing_momenta, 2)-1
    allocate(permutation(degree))
    permutation = [(order, order=1, degree)]
    do order = 1, size(amplitudes)
      call evaluate_one_order(outgoing_momenta, external_wavefunctions, &
           permutation, coupling, amplitudes(order))
      if (order < size(amplitudes)) call next_permutation(permutation)
    end do
  end subroutine evaluate_legacy_trace_orders

  subroutine evaluate_one_order(momenta, wavefunctions, permutation, &
       coupling, amplitude)
    real(dp), intent(in) :: momenta(0:, :)
    complex(dp), intent(in) :: wavefunctions(:, :)
    integer, intent(in) :: permutation(:)
    real(dp), intent(in) :: coupling
    complex(dp), intent(out) :: amplitude

    complex(dp), allocatable :: gluon(:, :, :), tensor(:, :, :)
    real(dp), allocatable :: interval(:, :, :)
    complex(dp) :: tensor_vertex(6), vector_vertex(4)
    integer :: degree, first, last, left_last, length, position

    degree = size(permutation)
    allocate(interval(0:3, degree, degree), gluon(4, degree, degree))
    allocate(tensor(6, degree, degree))
    interval = 0.0_dp
    gluon = cmplx(0.0_dp, 0.0_dp, dp)
    tensor = cmplx(0.0_dp, 0.0_dp, dp)
    do position = 1, degree
      interval(:, position, position) = momenta(:, permutation(position))
      gluon(:, position, position) = wavefunctions(:, permutation(position))
    end do

    do length = 2, degree
      do first = 1, degree-length+1
        last = first+length-1
        interval(:, first, last) = interval(:, first, last-1)+ &
             interval(:, last, last)
        do left_last = first, last-1
          call three_gluon(gluon(:, first, left_last), &
               interval(:, first, left_last), &
               gluon(:, left_last+1, last), &
               interval(:, left_last+1, last), vector_vertex)
          gluon(:, first, last) = gluon(:, first, last)+ &
               coupling*vector_vertex
          if (left_last > first) then
            call aux_tensor_gluon_to_gluon(tensor(:, first, left_last), &
                 gluon(:, left_last+1, last), vector_vertex)
            gluon(:, first, last) = gluon(:, first, last)+ &
                 coupling*vector_vertex
          end if
          if (left_last+1 < last) then
            call gluon_aux_tensor_to_gluon(gluon(:, first, left_last), &
                 tensor(:, left_last+1, last), vector_vertex)
            gluon(:, first, last) = gluon(:, first, last)+ &
                 coupling*vector_vertex
          end if
          if (length < degree) then
            call two_gluon_to_aux_tensor(gluon(:, first, left_last), &
                 gluon(:, left_last+1, last), tensor_vertex)
            tensor(:, first, last) = tensor(:, first, last)+ &
                 coupling*tensor_vertex
          end if
        end do
        if (length < degree) call massless_vector_propagator( &
             gluon(:, first, last), interval(:, first, last))
      end do
    end do
    amplitude = terminal_vector_contract(gluon(:, 1, degree), &
         wavefunctions(:, degree+1))
  end subroutine evaluate_one_order

  subroutine next_permutation(permutation)
    integer, intent(inout) :: permutation(:)

    integer :: left, right, temporary

    left = size(permutation)-1
    do while (permutation(left) >= permutation(left+1))
      left = left-1
    end do
    right = size(permutation)
    do while (permutation(right) < permutation(left))
      right = right-1
    end do
    temporary = permutation(left)
    permutation(left) = permutation(right)
    permutation(right) = temporary
    permutation(left+1:) = permutation(size(permutation):left+1:-1)
  end subroutine next_permutation

end module legacy_trace_order
