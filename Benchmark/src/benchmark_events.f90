module benchmark_events
  use, intrinsic :: iso_fortran_env, only: error_unit, real64
  implicit none
  private

  integer, parameter, public :: benchmark_dp = real64
  character(len=*), parameter :: event_magic = 'AMPLIGLUON_EVENT_V1'

  type, public :: benchmark_event_t
    integer :: final_gluons = 0
    integer :: number_of_helicities = 0
    real(benchmark_dp) :: strong_coupling = 1.0_benchmark_dp
    real(benchmark_dp), allocatable :: momenta(:, :)
    integer, allocatable :: helicities(:, :)
  contains
    procedure :: load => load_benchmark_event
  end type benchmark_event_t

  public :: benchmark_fail

contains

  subroutine load_benchmark_event(this, filename)
    class(benchmark_event_t), intent(inout) :: this
    character(len=*), intent(in) :: filename

    character(len=64) :: keyword
    integer, allocatable :: helicity_row(:)
    integer :: configuration, external_leg, input_unit, ios, total_gluons
    real(benchmark_dp) :: momentum_row(0:3)

    open(newunit=input_unit, file=filename, status='old', action='read', &
         form='formatted', iostat=ios)
    if (ios /= 0) call benchmark_fail('cannot open event file: '//trim(filename))

    read(input_unit, *, iostat=ios) keyword
    if (ios /= 0 .or. trim(keyword) /= event_magic) &
         call benchmark_fail('invalid event-file magic')
    read(input_unit, *, iostat=ios) keyword, this%final_gluons
    if (ios /= 0 .or. trim(keyword) /= 'FINAL_GLUONS' .or. &
         this%final_gluons < 2) &
         call benchmark_fail('invalid FINAL_GLUONS row')
    total_gluons = this%final_gluons+2
    read(input_unit, *, iostat=ios) keyword, this%strong_coupling
    if (ios /= 0 .or. trim(keyword) /= 'STRONG_COUPLING') &
         call benchmark_fail('invalid STRONG_COUPLING row')

    read(input_unit, *, iostat=ios) keyword
    if (ios /= 0 .or. trim(keyword) /= 'BEGIN_MOMENTA') &
         call benchmark_fail('missing BEGIN_MOMENTA row')
    if (allocated(this%momenta)) deallocate(this%momenta)
    allocate(this%momenta(0:3, total_gluons))
    do external_leg = 1, total_gluons
      read(input_unit, *, iostat=ios) momentum_row
      if (ios /= 0) call benchmark_fail('invalid momentum row')
      this%momenta(:, external_leg) = momentum_row
    end do
    read(input_unit, *, iostat=ios) keyword
    if (ios /= 0 .or. trim(keyword) /= 'END_MOMENTA') &
         call benchmark_fail('missing END_MOMENTA row')

    read(input_unit, *, iostat=ios) keyword, this%number_of_helicities
    if (ios /= 0 .or. trim(keyword) /= 'NHELICITIES' .or. &
         this%number_of_helicities < 1) &
         call benchmark_fail('invalid NHELICITIES row')
    read(input_unit, *, iostat=ios) keyword
    if (ios /= 0 .or. trim(keyword) /= 'BEGIN_HELICITIES') &
         call benchmark_fail('missing BEGIN_HELICITIES row')
    if (allocated(this%helicities)) deallocate(this%helicities)
    allocate(this%helicities(total_gluons, this%number_of_helicities))
    allocate(helicity_row(total_gluons))
    do configuration = 1, this%number_of_helicities
      read(input_unit, *, iostat=ios) helicity_row
      if (ios /= 0) call benchmark_fail('invalid helicity row')
      this%helicities(:, configuration) = helicity_row
    end do
    read(input_unit, *, iostat=ios) keyword
    if (ios /= 0 .or. trim(keyword) /= 'END_HELICITIES') &
         call benchmark_fail('missing END_HELICITIES row')
    close(input_unit)
  end subroutine load_benchmark_event

  subroutine benchmark_fail(message)
    character(len=*), intent(in) :: message

    write(error_unit, '(a)') 'Benchmark error: '//trim(message)
    error stop 1
  end subroutine benchmark_fail

end module benchmark_events
