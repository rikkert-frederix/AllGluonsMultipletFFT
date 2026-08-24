module rambo_on_diet
  ! Implementation of the 'RAMBO on diet' algorithm (Simon Plaetzer,
  ! arXiv:1308.2922 [hep-ph]). There are two versions of the momenta
  ! generation from the random numbers; the first is the
  ! implementation as is in the paper, the second uses instead the
  ! original 'RAMBO' algorithm to add the masses to the massless
  ! momenta (with the massless momenta generated according to the
  ! 'RAMBO on diet' algorithm). The first version is more efficient
  ! (weights closer to unity) when including masses. This is also the
  ! only version for which the inversion (random numbers from momenta)
  ! is implemented.
  use find_zero, only: find_zero_rambo_polynomial
  private
  real(kind=8),parameter :: pi=3.1415926535897932385d0
  public :: generate_mom_from_x_v1,generate_mom_from_x_v2, &
            reconstruct_x_from_mom_v1
contains
  subroutine generate_mom_from_x_v1(n,x,mass,Q_in,p,wgt)
    implicit none
    ! arguments
    integer,intent(in) :: n
    real(kind=8),intent(in) :: Q_in
    real(kind=8),dimension(n),intent(in) :: mass
    real(kind=8),dimension(3*n-4),intent(in) :: x
    real(kind=8),dimension(0:3,n),intent(out) :: p
    real(kind=8),intent(out) :: wgt
    ! local
    integer :: i
    real(kind=8),dimension(0:3,n) :: Q
    real(kind=8),dimension(n) :: M,K,u
    real(kind=8),dimension(1:3) :: v
    real(kind=8) :: q_tmp,costheta,sintheta,phi
    ! external
    real(kind=8) :: r
    Q(0,1)=Q_in
    Q(1:3,1)=0d0
    K(1)=Q_in-sum(mass(1:n))
    if (K(1).le.0d0) then
       write (*,*) 'ERROR: not enough CM energy to generate momenta:', &
            Q_in,'=< sum(',mass(1:n),')'
       stop 1
    endif
    K(n)=0d0
    do i=2,n-1
       r=x(i-1)
       call find_zero_rambo_polynomial(n-i,r,u(i))
       if (u(i).le.0d0 .or. u(i).ge.1d0) then
          write (*,*) 'ERROR: did not find a viable zero',i,u(i)
          stop 1
       endif
       K(i)=sqrt(u(i))*K(i-1)
    enddo
    do i=1,n
       M(i)=K(i)+sum(mass(i:n))
    enddo
    do i=2,n
       q_tmp=M(i-1)*rho_m(M(i-1),M(i),mass(i-1))/2d0
       costheta=2*x(n-5+2*i)-1
       sintheta=sqrt(1-costheta**2)
       phi=2d0*pi*x(n-4+2*i)
       p(1,i-1)=q_tmp*cos(phi)*sintheta
       p(2,i-1)=q_tmp*sin(phi)*sintheta
       p(3,i-1)=q_tmp*costheta
       p(0,i-1)=sqrt(q_tmp**2+mass(i-1)**2)
       v(1:3)=-Q(1:3,i-1)/Q(0,i-1)
       call boost(p(0,i-1),v,p(0,i-1))
       Q(1:3,i)=Q(1:3,i-1)-p(1:3,i-1)
       Q(0,i)=sqrt(M(i)**2+Q(1,i)**2+Q(2,i)**2+Q(3,i)**2)
    enddo
    p(0:3,n)=Q(0:3,n)
    wgt=(pi/2d0)**(n-1)*K(1)**(2*n-4)/(factorial(n-2)**2*(n-1))*(2d0*pi)**(4-3*n)
    do i=2,n-1
       wgt=wgt * (rho_m(M(i-1),M(i),mass(i-1))/rho_m(K(i-1),K(i),0d0)) * (M(i)/K(i))
    enddo
    wgt=wgt*rho_m(M(n-1),M(n),mass(n-1))
  end subroutine generate_mom_from_x_v1

  subroutine generate_mom_from_x_v2(n,x,mass,Q_in,p,wgt)
    implicit none
    integer,intent(in) :: n
    real(kind=8),intent(in) :: Q_in
    real(kind=8),dimension(n),intent(in) :: mass
    real(kind=8),dimension(3*n-4),intent(in) :: x
    real(kind=8),dimension(0:3,n),intent(out) :: p
    real(kind=8),intent(out) :: wgt
    real(kind=8),dimension(0:3,n) :: pm
    call generate_massless_mom_from_x(n,x,Q_in,p,wgt)
    call add_masses_to_massless(n,Q_in,p,mass,pm,wgt)
    p=pm
  end subroutine generate_mom_from_x_v2

  
  subroutine generate_massless_mom_from_x(n,x,Q_in,p,wgt)
    implicit none
    ! arguments
    integer,intent(in) :: n
    real(kind=8),intent(in) :: Q_in
    real(kind=8),dimension(3*n-4),intent(in) :: x
    real(kind=8),dimension(0:3,n),intent(out) :: p
    real(kind=8),intent(out) :: wgt
    ! local
    real(kind=8),dimension(n) :: mass
    mass=0d0
    call generate_mom_from_x_v1(n,x,mass,Q_in,p,wgt)
  end subroutine generate_massless_mom_from_x

  subroutine add_masses_to_massless(n,Q_in,p,mass,pm,wgt)
    implicit none
    ! argument
    integer,intent(in) :: n
    real(kind=8),intent(in) :: Q_in
    real(kind=8),dimension(n),intent(in) :: mass
    real(kind=8),dimension(0:3,n),intent(in) :: p
    real(kind=8),dimension(0:3,n),intent(out) :: pm
    real(kind=8),intent(inout) :: wgt
    ! local
    integer :: iter,i
    real(kind=8) :: mass_tot,xi,f0,g0,xi2,accu,wgt2,wgt3
    real(kind=8),dimension(n) :: mass2,p2,E,v
    ! constant
    integer,parameter :: itmax=6
    real(kind=8),parameter :: acc=1d-14
    mass_tot=sum(mass(1:n))
    if (abs(mass_tot).le.tiny(1d0)) then
       pm=p
       return
    endif
    do i=1,n
       mass2(i)=mass(i)**2
       p2(i)=p(0,i)**2
    enddo
    iter=0
    xi=sqrt(1d0-(mass_tot/Q_in)**2)
    accu=Q_in*acc
    do
       f0=-Q_in
       g0=0d0
       xi2=xi**2
       do i=1,n
          E(i)=sqrt(mass2(i)+xi2*p2(i))
          f0=f0+E(i)
          g0=g0+p2(i)/E(i)
       enddo
       if (abs(f0).le.accu) exit
       iter=iter+1
       if (iter.eq.itmax) then
          write (*,*) 'lack of precision in masses',f0,accu
          stop 1
       endif
       xi=xi-f0/(xi*g0)
    enddo
    do i=1,n
       pm(1:3,i)=xi*p(1:3,i)
       pm(0,i)=E(i)
    enddo
    wgt2=1d0
    wgt3=0d0
    do i=1,n
       v(i)=xi*p(0,i)
       wgt2=wgt2*v(i)/E(i)
       wgt3=wgt3+v(i)**2/E(i)
    enddo
    wgt=wgt*xi**(2*n-3)*wgt2/wgt3*Q_in
  end subroutine add_masses_to_massless
  

  subroutine reconstruct_x_from_mom_v1(n,p,mass,Q_in,x)
    implicit none
    ! arguments
    integer,intent(in) :: n
    real(kind=8),intent(in) :: Q_in
    real(kind=8),dimension(n),intent(in) :: mass
    real(kind=8),dimension(3*n-4),intent(out) :: x
    real(kind=8),dimension(0:3,n),intent(in) :: p
    ! local
    integer :: i
    real(kind=8),dimension(0:3,n) :: Q
    real(kind=8),dimension(n) :: M,K
    real(kind=8),dimension(1:3) :: v
    real(kind=8),dimension(0:3) :: p_tot
    real(kind=8) :: phi,u
    M(1)=Q_in
    Q(0:3,n)=p(0:3,n)
    p_tot(0:3)=0d0
    M(n)=0d0
    do i=n,2,-1
       p_tot(0:3)=p_tot(0:3)+p(0:3,i)
       if (i.ne.n) m(i)=sqrt(p_tot(0)**2-p_tot(1)**2-p_tot(2)**2-p_tot(3)**2)
    enddo
    do i=1,n
       K(i)=M(i)-sum(mass(i:n))
    enddo
    do i=n,2,-1
       u=(K(i)/K(i-1))**2
       if (i.ne.n) x(i-1)=(n+1-i)*u**(n-i)-(n-i)*u**(n+1-i)
       Q(0:3,i-1)=Q(0:3,i)+p(0:3,i-1)
       v(1:3)=Q(1:3,i-1)/Q(0,i-1)
       call boost(p(0,i-1),v,p(0,i-1))
       x(n-5+2*i)=(p(3,i-1)/sqrt(p(1,i-1)**2+p(2,i-1)**2+p(3,i-1)**2)+1d0)/2d0
       phi=atan2(p(2,i-1),p(1,i-1))
       if (phi.lt.0d0) phi=phi+2*pi
       x(n-4+2*i)=phi/(2*pi)
    enddo
  end subroutine reconstruct_x_from_mom_v1

  subroutine boost(p,v,q)
    implicit none
    real(kind=8),dimension(0:3) :: p,q
    real(kind=8),dimension(1:3) :: v,v_norm
    real(kind=8) :: beta,gamma,chybst,shybst,chybstmo,pv,en
    beta=sqrt(v(1)**2+v(2)**2+v(3)**2)
    if (beta.le.1d-8) then
       q(0:3)=p(0:3)
       return
    endif
    v_norm(1:3)=v(1:3)/beta
    gamma=1d0/sqrt(1d0-beta**2)
    chybst=gamma
    shybst=gamma*beta
    chybstmo=gamma-1d0
    en=p(0)
    pv=p(1)*v_norm(1)+p(2)*v_norm(2)+p(3)*v_norm(3)
    q(0)=en*chybst-pv*shybst
    q(1:3)=p(1:3)+v_norm(1:3)*(pv*chybstmo-en*shybst)
  end subroutine boost
  
  real(kind=8) function rho_m(m1,m2,m3)
    implicit none
    real(kind=8) :: m1,m2,m3,tmp1,tmp2
    tmp1=(m1-(m2+m3))*(m1+(m2+m3))
    tmp2=(m1-(m2-m3))*(m1+(m2-m3))
    rho_m=sqrt(tmp1*tmp2)/m1**2
  end function rho_m

  real(kind=8) function factorial(ifact)
    ! computes the factorial of 'ifact'
    implicit none
    integer, value :: ifact
    factorial=1d0
    do while (ifact.gt.1)
       factorial=factorial*real(ifact,kind=8)
       ifact=ifact-1
    enddo
  end function factorial

end module rambo_on_diet
