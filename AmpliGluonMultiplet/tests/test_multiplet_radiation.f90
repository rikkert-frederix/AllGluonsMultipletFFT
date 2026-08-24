program test_multiplet_radiation
  use ampligluon_multiplet_kinds, only: dp, fail
  use multiplet_paths, only: path_catalog_t
  use multiplet_radiation, only: multiplet_radiation_system_t, &
       multiplet_radiation_workspace_t
  use recoupling_plan, only: recoupling_system_t
  use wigner_table, only: wigner_table_t
  implicit none

  integer, parameter :: maximum_length = 7
  real(dp), parameter :: structural_tolerance = 2.0e-12_dp

  type(path_catalog_t) :: catalog
  type(multiplet_radiation_system_t) :: radiation
  type(multiplet_radiation_workspace_t) :: workspace
  type(recoupling_system_t) :: recoupling
  type(wigner_table_t) :: table
  character(len=1024) :: table_file
  complex(dp), allocatable :: emitted(:), source(:), summed(:)
  real(dp) :: residual
  integer :: emitter, new_length, old_length, source_path

  if (command_argument_count() /= 1) &
       call fail('usage: test_multiplet_radiation WIGNER_TABLE')
  call get_command_argument(1, table_file)
  call table%load(trim(table_file))
  call catalog%build(table, maximum_length)
  call recoupling%build(table, catalog)
  call radiation%build(catalog)

  do old_length = 1, maximum_length-1
    new_length = old_length+1
    call require(radiation%source_dimension(old_length) == &
                 catalog%spaces(old_length)%number_of_paths, &
                 'wrong radiation source dimension')
    call require(radiation%target_dimension(old_length) == &
                 catalog%spaces(new_length)%number_of_paths, &
                 'wrong radiation target dimension')
    call require(radiation%number_of_emitters(old_length) == old_length+1, &
                 'wrong number of radiation emitters')

    allocate(source(radiation%source_dimension(old_length)))
    allocate(emitted(radiation%target_dimension(old_length)))
    allocate(summed(radiation%target_dimension(old_length)))

    ! Cover every column of every radiation operator.  Colour conservation is
    ! a tensor identity, so checking each old-basis unit vector is a complete
    ! structural check rather than a random numerical sample.
    residual = 0.0_dp
    do source_path = 1, size(source)
      source = (0.0_dp, 0.0_dp)
      source(source_path) = cmplx(1.0_dp, -0.375_dp, dp)
      summed = (0.0_dp, 0.0_dp)
      do emitter = 1, radiation%number_of_emitters(old_length)
        call radiation%add(recoupling, old_length, emitter, source, &
                           cmplx(1.0_dp, 0.0_dp, dp), summed, workspace)
      end do
      residual = max(residual, maxval(abs(summed)))
    end do
    call require(residual <= structural_tolerance, &
                 'sum of colour-charge insertions is nonzero')

    ! Each construction is an isometry: adjacent swaps are orthogonal and
    ! appending the normalized multiplicity-one vertex is one-to-one.
    do source_path = 1, size(source)
      source(source_path) = cmplx(cos(real(source_path, dp)), &
                                  sin(real(2*source_path, dp)), dp)
    end do
    do emitter = 1, radiation%number_of_emitters(old_length)
      emitted = cmplx(-7.0_dp, 3.0_dp, dp)
      call radiation%apply(recoupling, old_length, emitter, source, &
                           emitted, workspace)
      call require(abs(squared_norm(emitted)-squared_norm(source)) <= &
                   structural_tolerance*max(1.0_dp, squared_norm(source)), &
                   'radiation operator is not isometric')
    end do

    deallocate(source, emitted, summed)
  end do

  write(*, '(a)') 'multiplet radiation regression: PASS'

contains

  real(dp) function squared_norm(vector) result(value)
    complex(dp), intent(in) :: vector(:)

    value = sum(real(vector*conjg(vector), dp))
  end function squared_norm

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) call fail(message)
  end subroutine require

end program test_multiplet_radiation
