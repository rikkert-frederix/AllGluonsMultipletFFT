module find_zero

  ! Finds the zero of a function between x=0d0 and x=1d0 for the
  ! function f(x) with deriviative fprime(x). The function should pass
  ! once through the zero in the interval [0,1].  Uses the Bisection
  ! method to find a value close to the zero, and then switches to
  ! Newton's method to zoom in on the correct result.

contains

  subroutine find_zero_rambo_polynomial(power,r,c)
    ! Solve (power+1)*c**power-power*c**(power+1)=r on (0,1).
    ! This specialized, bracketed form avoids the executable-stack
    ! trampoline needed when an internal procedure is passed as a callback.
    implicit none

    integer, intent(in) :: power
    real(kind=8), intent(in) :: r
    real(kind=8), intent(out) :: c

    integer :: iteration
    real(kind=8) :: lower, upper, value

    if (power.lt.1) then
       write (*,*) 'Error: RAMBO polynomial power must be positive',power
       stop 1
    endif
    if (r.le.0d0 .or. r.ge.1d0) then
       write (*,*) 'Error: RAMBO random number must lie in (0,1)',r
       stop 1
    endif

    lower=0d0
    upper=1d0
    do iteration=1,60
       c=(lower+upper)/2d0
       value=(power+1)*c**power-power*c**(power+1)-r
       if (value.lt.0d0) then
          lower=c
       else
          upper=c
       endif
    enddo
    c=(lower+upper)/2d0
  end subroutine find_zero_rambo_polynomial

  subroutine find_zero_pol(f,fprime,c)
  
  implicit none

  ! Tolerance for the zero; the code continues to find the zero x_zero
  ! until f(x_zero) < tol.
  real(kind=8),parameter :: tol=1d-12

  ! always perform at least 'min_bisect' bisection iteration (except
  ! if we get real lucky and find a zero within those
  ! iterations). Perform more iterations if we keep getting that the
  ! zero is in the interval next to the boundary, up to a maximum of
  ! 'max_bisect'. In case the zero has not been found, use up to
  ! 'max_newton' iterations of Newton's method.
  integer,parameter :: min_bisect=6
  integer,parameter :: max_bisect=40
  integer,parameter :: max_newton=40

  ! local variables:
  real(kind=8) :: a,b,fa,fb,tmp,res,der,c
  logical :: found
  integer :: k,iden_size
  
! The function (see below)
  real(kind=8) :: f,fprime

  ! found zero between these two boundaries:
  a=0d0
  b=1d0

  fa=f(a)
  fb=f(b)
  if ( (fa.lt.0d0) .and. (fb.gt.0d0) ) then
     ! Great. Have an increasing function passing through zero
     continue
  elseif ( (fa.gt.0d0) .and. (fb.lt.0d0) ) then
     ! Still great, but need to switch the boundaries:
     tmp=a
     a=b
     b=tmp
  else
     write (*,*) 'Error: Cannot found zero because function might not pass through zero',fa,fb
     stop 1
  endif

  ! Start by using the bisection method to get a value close to the
  ! true value:
  found=.false.
  iden_size=0
  do k=3,max_bisect ! start counter at '3', since first two are the boundary values
     c=(a+b)/2
     res=f(c)
!!$     write (*,*) 'bisect',k,c,res
     if(abs(res).lt.tol) then
        ! in the rare case that we find the zero already here, quit
        ! the loop before it finishes
        found=.true.
        exit
     endif
     if(res.LT.0)THEN
        a=c
        iden_size=iden_size+1
     else
        b=c
        iden_size=iden_size-1
     endif

     if (k.ne.abs(iden_size)+2 .and. k.ge.min_bisect) then
        ! Exit the loop if we are no longer at the intervals that is
        ! next to the boundary and we have done a minimum number of
        ! iterations
        exit
     endif
  enddo

  ! If not yet found (which should be effectively always), use
  ! Newton's method to refine the results:
  if (.not. found) then
     der=fprime(c)
     c=c-res/der
     do k=1,max_newton
        res=f(c)
!!$        write (*,*) 'newton',k,c,res
        if (abs(res).lt.tol) then
           ! found the zero; exit the loop
           found=.true.
           exit
        endif
        der=fprime(c)
        c=c-res/der
     enddo
  endif

  if (found) then
!!$     write (*,*) 'zero is:',c,res
  else
     write (*,*) 'zero not found'
     stop 1
  endif

end subroutine find_zero_pol

  
end module find_zero
