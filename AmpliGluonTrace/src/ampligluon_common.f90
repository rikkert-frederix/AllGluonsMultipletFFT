module ampligluon_common
  use, intrinsic :: iso_fortran_env, only: error_unit, int16, int64, real64
  implicit none
  private

  integer, parameter, public :: dp = real64
  integer, parameter, public :: i16 = int16
  integer, parameter, public :: i64 = int64

  public :: fail

contains

  subroutine fail(message)
    character(len=*), intent(in) :: message

    write(error_unit, '(a)') 'All-gluon amplitude error: '//trim(message)
    error stop 1
  end subroutine fail

end module ampligluon_common
