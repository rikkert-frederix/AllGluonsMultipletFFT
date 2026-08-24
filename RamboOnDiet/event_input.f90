module event_input
  use ampligluon_multiplet_kinds, only: dp, fail
  implicit none
  private

  character(len=*), parameter :: event_magic = 'AMPLIGLUON_EVENT_V1'

  type, public :: event_input_t
    integer :: final_gluons = 0
    integer :: number_of_helicities = 0
    real(dp) :: strong_coupling = 1.0_dp
    real(dp), allocatable :: momenta(:, :)
    integer, allocatable :: helicities(:, :)
  contains
    procedure :: load => load_event
  end type event_input_t

contains

  subroutine load_event(this, filename)
    class(event_input_t), intent(inout) :: this
    character(len=*), intent(in) :: filename

    character(len=64) :: keyword
    integer, allocatable :: helicity_row(:)
    integer :: configuration, external_leg, input_unit, ios, total_gluons
    real(dp) :: momentum_row(0:3)

    open(newunit=input_unit, file=filename, status='old', action='read', &
         form='formatted', iostat=ios)
    if (ios /= 0) call fail('cannot open event file: '//trim(filename))

    read(input_unit, *, iostat=ios) keyword
    if (ios /= 0 .or. trim(keyword) /= event_magic) &
         call fail('invalid event-file magic')
    read(input_unit, *, iostat=ios) keyword, this%final_gluons
    if (ios /= 0 .or. trim(keyword) /= 'FINAL_GLUONS' .or. &
        this%final_gluons < 2) call fail('invalid FINAL_GLUONS row')
    total_gluons = this%final_gluons+2
    read(input_unit, *, iostat=ios) keyword, this%strong_coupling
    if (ios /= 0 .or. trim(keyword) /= 'STRONG_COUPLING') &
         call fail('invalid STRONG_COUPLING row')

    read(input_unit, *, iostat=ios) keyword
    if (ios /= 0 .or. trim(keyword) /= 'BEGIN_MOMENTA') &
         call fail('missing BEGIN_MOMENTA row')
    if (allocated(this%momenta)) deallocate(this%momenta)
    allocate(this%momenta(0:3, total_gluons))
    do external_leg = 1, total_gluons
      read(input_unit, *, iostat=ios) momentum_row
      if (ios /= 0) call fail('invalid momentum row')
      this%momenta(:, external_leg) = momentum_row
    end do
    read(input_unit, *, iostat=ios) keyword
    if (ios /= 0 .or. trim(keyword) /= 'END_MOMENTA') &
         call fail('missing END_MOMENTA row')

    read(input_unit, *, iostat=ios) keyword, this%number_of_helicities
    if (ios /= 0 .or. trim(keyword) /= 'NHELICITIES' .or. &
        this%number_of_helicities < 1) call fail('invalid NHELICITIES row')
    read(input_unit, *, iostat=ios) keyword
    if (ios /= 0 .or. trim(keyword) /= 'BEGIN_HELICITIES') &
         call fail('missing BEGIN_HELICITIES row')
    if (allocated(this%helicities)) deallocate(this%helicities)
    allocate(this%helicities(total_gluons, this%number_of_helicities))
    allocate(helicity_row(total_gluons))
    do configuration = 1, this%number_of_helicities
      read(input_unit, *, iostat=ios) helicity_row
      if (ios /= 0) call fail('invalid helicity row')
      this%helicities(:, configuration) = helicity_row
    end do
    read(input_unit, *, iostat=ios) keyword
    if (ios /= 0 .or. trim(keyword) /= 'END_HELICITIES') &
         call fail('missing END_HELICITIES row')
    close(input_unit)
  end subroutine load_event

end module event_input
