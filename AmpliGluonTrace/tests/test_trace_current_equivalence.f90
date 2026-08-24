program test_trace_current_equivalence
  use ampligluon_common, only: dp, fail
  use gluon_kinematics, only: external_massless_vector
  use legacy_trace_order, only: evaluate_legacy_trace_orders
  use trace_current_dag, only: trace_current_dag_t
  use trace_mhv, only: trace_mhv_t
  use trace_order_recursion, only: trace_order_recursion_t
  implicit none

  integer :: total_gluons

  do total_gluons = 4, 9
    call compare_all_orders(total_gluons)
  end do
  write(*, '(a)') 'trace current equivalence regression: PASS'

contains

  subroutine compare_all_orders(total_gluons)
    integer, intent(in) :: total_gluons

    type(trace_current_dag_t) :: dag, fixed_dag
    type(trace_mhv_t) :: mhv
    type(trace_order_recursion_t) :: order_recursion
    complex(dp), allocatable :: external(:, :), fallback(:), fixed(:), joint(:)
    complex(dp), allocatable :: legacy(:), mhv_values(:)
    real(dp), allocatable :: outgoing(:, :), physical(:, :)
    real(dp) :: difference, scale
    real(dp), parameter :: couplings(3) = [1.0_dp, 0.73_dp, -0.5_dp]
    integer :: gluon, helicities(total_gluons), order

    allocate(physical(0:3, total_gluons), outgoing(0:3, total_gluons))
    allocate(external(4, total_gluons))
    call fill_momenta(total_gluons, physical)
    outgoing(:, 1:2) = -physical(:, 1:2)
    outgoing(:, 3:) = physical(:, 3:)
    helicities = [(merge(-1, 1, modulo(gluon, 2) == 1), &
         gluon=1, total_gluons)]
    do gluon = 1, total_gluons
      call external_massless_vector(outgoing(:, gluon), helicities(gluon), &
           external(:, gluon))
    end do

    call dag%initialize(total_gluons-1)
    allocate(fallback(dag%number_of_orders()), joint(dag%number_of_orders()))
    allocate(legacy(dag%number_of_orders()))
    allocate(mhv_values(dag%number_of_orders()))
    call mhv%initialize(total_gluons, dag%number_of_orders())
    call fixed_dag%initialize(total_gluons-1, fixed_first=.true.)
    allocate(fixed(fixed_dag%number_of_orders()))
    call dag%evaluate(outgoing(:, :total_gluons-1), &
         external(:, :total_gluons-1), external(:, total_gluons), &
         0.73_dp, joint)
    call fixed_dag%evaluate(outgoing(:, :total_gluons-1), &
         external(:, :total_gluons-1), external(:, total_gluons), &
         0.73_dp, fixed)
    call evaluate_legacy_trace_orders(outgoing, external, 0.73_dp, legacy)
    call order_recursion%initialize(total_gluons-1, dag%number_of_orders())
    call order_recursion%evaluate(outgoing(:, :total_gluons-1), &
         external(:, :total_gluons-1), external(:, total_gluons), &
         0.73_dp, fallback)
    do order = 1, size(joint)
      difference = abs(joint(order)-legacy(order))
      scale = max(abs(joint(order)), abs(legacy(order)), tiny(1.0_dp))
      if (difference > 2.0e-11_dp*scale) then
        write(*, '(a,2i0,3(1x,es24.16))') &
             'ordered-current mismatch at N/order ', total_gluons, order, &
             joint(order), legacy(order), difference/scale
        call fail('shared and legacy trace currents disagree')
      end if
      difference = abs(fallback(order)-legacy(order))
      scale = max(abs(fallback(order)), abs(legacy(order)), tiny(1.0_dp))
      if (difference > 2.0e-11_dp*scale) &
           call fail('fallback and legacy trace currents disagree')
    end do
    do order = 1, size(fixed)
      difference = abs(fixed(order)-legacy(order))
      scale = max(abs(fixed(order)), abs(legacy(order)), tiny(1.0_dp))
      if (difference > 2.0e-11_dp*scale) then
        write(*, '(a,2i0,3(1x,es24.16))') &
             'fixed-first current mismatch at N/order ', total_gluons, order, &
             fixed(order), legacy(order), difference/scale
        call fail('fixed-first and legacy trace currents disagree')
      end if
    end do

    do order = 1, size(couplings)
      helicities = -1
      call compare_mhv_configuration(total_gluons, outgoing, helicities, &
           external, dag, mhv, joint, mhv_values, couplings(order))
      helicities = 1
      call compare_mhv_configuration(total_gluons, outgoing, helicities, &
           external, dag, mhv, joint, mhv_values, couplings(order))
      if (total_gluons >= 5) then
        helicities = 1
        helicities(1:2) = -1
        helicities(3:4) = -1
        call compare_mhv_configuration(total_gluons, outgoing, helicities, &
             external, dag, mhv, joint, mhv_values, couplings(order))
        helicities = -1
        helicities(1:2) = 1
        helicities(3:4) = 1
        call compare_mhv_configuration(total_gluons, outgoing, helicities, &
             external, dag, mhv, joint, mhv_values, couplings(order))
      end if
    end do
  end subroutine compare_all_orders

  subroutine compare_mhv_configuration(total_gluons, outgoing, helicities, &
       external, dag, mhv, dag_values, mhv_values, coupling)
    integer, intent(in) :: total_gluons, helicities(total_gluons)
    real(dp), intent(in) :: outgoing(0:3, total_gluons)
    complex(dp), intent(inout) :: external(4, total_gluons)
    type(trace_current_dag_t), intent(inout) :: dag
    type(trace_mhv_t), intent(inout) :: mhv
    complex(dp), intent(inout) :: dag_values(:), mhv_values(:)
    real(dp), intent(in) :: coupling

    real(dp) :: difference, scale
    integer :: gluon, order

    do gluon = 1, total_gluons
      call external_massless_vector(outgoing(:, gluon), helicities(gluon), &
           external(:, gluon))
    end do
    call dag%evaluate(outgoing(:, :total_gluons-1), &
         external(:, :total_gluons-1), external(:, total_gluons), &
         coupling, dag_values)
    call mhv%evaluate(outgoing, external, helicities, coupling, mhv_values)
    do order = 1, size(dag_values)
      difference = abs(mhv_values(order)-dag_values(order))
      scale = max(abs(mhv_values(order)), abs(dag_values(order)), tiny(1.0_dp))
      if (difference > 5.0e-11_dp*scale) then
        write(*, '(a,2i0,3(1x,es24.16))') &
             'MHV trace-order mismatch at N/order ', total_gluons, order, &
             mhv_values(order), dag_values(order), difference/scale
        call fail('Parke-Taylor and shared trace currents disagree')
      end if
    end do
  end subroutine compare_mhv_configuration

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
    end select
  end subroutine fill_momenta

end program test_trace_current_equivalence
