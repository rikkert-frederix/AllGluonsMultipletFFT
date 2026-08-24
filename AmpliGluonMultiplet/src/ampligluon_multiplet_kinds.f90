module ampligluon_multiplet_kinds
  use, intrinsic :: iso_fortran_env, only: error_unit, int64, real64
  implicit none
  private

  integer, parameter, public :: dp = real64
  integer, parameter, public :: i8 = int64
  integer, parameter, public :: stderr = error_unit
  complex(dp), parameter, public :: imaginary_unit = (0.0_dp, 1.0_dp)

  public :: fail

contains

  subroutine fail(message)
    character(len=*), intent(in) :: message

    write(stderr, '(a)') 'AmpliGluonMultiplet: '//trim(message)
    error stop 1
  end subroutine fail

end module ampligluon_multiplet_kinds
