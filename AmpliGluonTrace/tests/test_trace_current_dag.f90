program test_trace_current_dag
  use ampligluon_common, only: fail
  use trace_current_dag, only: trace_current_dag_t
  implicit none

  integer, parameter :: number_of_cases = 6
  integer, parameter :: degrees(number_of_cases) = [3, 4, 5, 6, 7, 8]
  integer, parameter :: expected_orders(number_of_cases) = &
       [6, 24, 120, 720, 5040, 40320]
  integer, parameter :: expected_word_slots(number_of_cases) = &
       [9, 34, 165, 981, 6853, 54804]
  integer, parameter :: expected_interactions(number_of_cases) = &
       [6, 33, 205, 1470, 11991, 109606]
  integer, parameter :: expected_attachments(number_of_cases) = &
       [9, 66, 490, 3915, 34251, 328804]
  integer, parameter :: expected_fixed_orders(number_of_cases) = &
       [2, 6, 24, 120, 720, 5040]
  integer, parameter :: expected_fixed_word_slots(number_of_cases) = &
       [8, 25, 99, 491, 2938, 20553]
  integer, parameter :: expected_fixed_interactions(number_of_cases) = &
       [6, 30, 163, 1020, 7341, 59941]
  integer, parameter :: expected_fixed_attachments(number_of_cases) = &
       [7, 42, 262, 1795, 13701, 116452]

  type(trace_current_dag_t) :: dag
  integer :: test_case

  do test_case = 1, number_of_cases
    call dag%initialize(degrees(test_case))
    call require(dag%number_of_orders() == expected_orders(test_case), &
         'wrong trace-current order count')
    call require(dag%number_of_word_slots() == expected_word_slots(test_case), &
         'wrong canonical trace-current word-slot count')
    call require(dag%number_of_interactions() == &
         expected_interactions(test_case), &
         'wrong shared trace interaction count')
    call require(dag%number_of_attachments() == &
         expected_attachments(test_case), &
         'wrong signed trace-current attachment count')

    call dag%initialize(degrees(test_case), fixed_first=.true.)
    call require(dag%number_of_orders() == expected_fixed_orders(test_case), &
         'wrong fixed-first trace-current order count')
    call require(dag%number_of_word_slots() == &
         expected_fixed_word_slots(test_case), &
         'wrong fixed-first trace-current word-slot count')
    call require(dag%number_of_interactions() == &
         expected_fixed_interactions(test_case), &
         'wrong fixed-first trace-current interaction count')
    call require(dag%number_of_attachments() == &
         expected_fixed_attachments(test_case), &
         'wrong fixed-first trace-current attachment count')
  end do

  write(*, '(a)') 'trace current DAG regression: PASS'

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) call fail(message)
  end subroutine require

end program test_trace_current_dag
