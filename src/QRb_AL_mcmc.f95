! Written by Dries F. Benoit 
! Faculty of economics and business administration
! Ghent University - BELGIUM
! LAST UPDATE (dd/mm/yy): 02/12/12

! MCMC sampler for binary quantile regression with 
! adaptive lasso variable selection.
! This code is and  

! Input arguments:
!    - n          : number of units of analysis
!    - k          : number of independent variables
!    - r          : number of MCMC iterations
!    - keep       : thinning parameter of MCMC draws
!    - y          : dependent variable
!    - p          : quantile of interest
!    - x          : matrix of regressors (1:n, 1:nvar)
!    - c          : shape parameter of gamma prior on etasquared 
!    - d          : rate parameter of gamma prior on etasquared 

! Output arguments:
!    - betadraw	  : the matrix of regression parameter estimates


subroutine QRb_AL_mcmc (n, k, r, keep, y, p, x, c, d, betadraw)

implicit none

! Precision statement:
integer, parameter :: dp = kind(1.0d0)

! Input arguments:
integer, intent(in) :: n, k, r, keep
integer, intent(in), dimension(n) :: y
real(dp), intent(in) :: p, c, d 
real(dp), intent(in), dimension(n,k) :: x

! Output arguments:
real(dp), intent(out), dimension(r/keep,k) :: betadraw

! Internal arguments:
integer :: i1, i2
real(dp) :: theta, omegasq, sigma, sigsq, shape, scale
real(dp), dimension(k) :: s, beta, etasq, lambda, beta0
real(dp), dimension(n) :: nu, mu, ytilde, ystar, Xbeta


! -- USEFUL QUANTITIES 
theta = (1.0_dp - 2.0_dp*p)/(p*(1.0_dp - p))
omegasq = 2.0_dp/(p*(1.0_dp - p))

! -- SET STARTING VALUES
ystar = 1.0_dp
nu = 1.0_dp
s = 1.0_dp
beta = 1.0_dp
sigma = 1.0_dp
etasq = 1.0_dp


! -- START OF MCMC CHAIN
do i1 = 1,r

  ! Simulate new values for ystar 
  !ooooooooooooooooooooooooooooooo
    Xbeta = matmul(X,beta)
    do i2 = 1,n
      if (y(i2) == 0) then
        call rtnorm_geweke(0.0_dp, .false., Xbeta(i2)+theta*nu(i2), sqrt(omegasq)&
                           *sqrt(nu(i2)), ystar(i2))
      else
        call rtnorm_geweke(0.0_dp, .true., Xbeta(i2)+theta*nu(i2), sqrt(omegasq)&
                           *sqrt(nu(i2)), ystar(i2))
      endif
    enddo

  ! Simulate new values for nu 
  !oooooooooooooooooooooooooooo
    lambda(1) = theta**2/(sigma*omegasq) + 2.0_dp/sigma
    mu = sqrt(lambda(1)/((ystar-matmul(x,beta))**2/(sigma*omegasq)))
    do i2 = 1,n
        call rinvgaus(mu(i2),lambda(1),nu(i2))
    enddo
    nu = 1.0_dp/nu

  ! Simulate new values for s 
  !ooooooooooooooooooooooooooo
    lambda = etasq
    mu(1:k) = sqrt(etasq/beta**2)
    do i2 = 1,k
        call rinvgaus(mu(i2),lambda(i2),s(i2))
    enddo 
    s = 1.0_dp/s

  ! Simulate new values for beta 
  !oooooooooooooooooooooooooooooo
    do i2 = 1,k
        beta0 = beta
        beta0(i2) = 0.0_dp
        ytilde = ystar - theta*nu - matmul(x,beta0)
        sigsq = 1.0_dp/(sum((x(:,i2)**2)/nu)/(sigma*omegasq) + 1.0_dp/s(i2))
        mu(1) = sigsq/(sigma*omegasq)*sum(ytilde*x(:,i2)/nu)
        call rnorm(beta(i2))
        beta(i2) = mu(1) + beta(i2)*sqrt(sigsq)
    enddo

  ! Simulate new values for etasq 
  !ooooooooooooooooooooooooooooooo
    shape = 1.0_dp + c
    do i2 = 1,k
        scale = 1.0_dp/(s(i2)/2.0_dp + d)
        if (scale < 1.0d-2) scale = 1.0d-2 !numerical stability
        if (scale > 1.0d2) scale = 1.0d2 !numerical stability
        call rgamma(shape,scale,etasq(i2))
    enddo

  ! Save current draw 
  !ooooooooooooooooooo
    if (mod(i1, keep) == 0) then
      betadraw((i1/keep),1:k) = beta
    endif

  ! Print information to console 
  !oooooooooooooooooooooooooooooo
    if (mod(i1, 500) == 0) then
      !write(*,*) 'Current iteration :', i1
      call intpr('Current iteration :', -1, i1, 1)
    endif

enddo

!===========================================================================================

contains

!===========================================================================================


! This code generates one draw from the standard normal 
! distribution. Note that more efficient code is possible
! when more than one normal draw is required.
! This code is based on the Box-Muller method.

! Output arguments:
!	- fn_val	: random draw from N(0,1) distribution

subroutine rnorm(fn_val)

implicit none

! Precision statement:
integer, parameter :: dp = kind(1.0d0)

! Output arguments:
real(dp), intent(out) :: fn_val

! Internal arguments:
real(dp) :: pi
real(dp), dimension(1:2) :: u

pi = 3.14159265358979323846_dp

call random_number(u)

fn_val = sqrt(-2*log(u(1))) * cos(2*pi*u(2))

end subroutine rnorm


!===========================================================================================


! Returns one draw from the truncated normal distribution

! Algorithm based on:
! Geweke, J. (1991). Efficient Simulation From the Multivariate Normal 
! and Student t-Distributions Subject to Linear Constraints, in Computer 
! Sciences and Statistics Proceedings of the 23d Symposium on the 
! Interface, pp. 571-578.

! This subroutine makes use of the subroutines:
!	- rnorm		: Box-Muller method for random normal draws

! Input arguments:
! a             -	trucation point
! lb	        -	logical:        if .TRUE. then trucation (a,+Inf)
!		        		if .FALSE. then truncation (-Inf,a)
! mu	        -	mean of trunc normal
! sigma         -	sd of trunc normal
! fn_val        -	random draw from trunc normal


subroutine rtnorm_geweke(a, lb, mu, sigma, fn_val)
  
implicit none

! Precision statement:
integer, parameter :: dp = kind(1.0d0)

! Input arguments:
logical, intent(in) :: lb
real(dp), intent(in) :: a, mu, sigma

! Output arguments:
real(dp), intent(out) :: fn_val

! Internal arguments:
real(dp) :: z, phi_z, az, c
real(dp), dimension(1:2) :: u


! Rescale truncation point
az=(a-mu)/sigma

if (lb) then
    c=az
  else 
    c=-az
endif

if (c<.45_dp) then

  ! normal rejection sampling
  do
    call rnorm(u(1))
    
    if (u(1)>c) exit
  end do
  z=u(1)

else

  ! exponential rejection sampling
  do
    ! Create exponential random variate z
    ! from uniform random variate u(1)
    call random_number(u)
    z = -log(u(1))/c

    phi_z = exp(-.5_dp * z**2_dp) !see Geweke
    if (u(2)<phi_z) exit
  end do
  z=z+c

end if

if (lb) then 
  fn_val = mu + sigma*z
else
  fn_val = mu - sigma*z
end if

end subroutine rtnorm_geweke


!===========================================================================================


! This code generates one random draw from the inverse Gaussian distribution.
! The algorithm is based on: Michael, Schucany & Haas (1976), Generating
! random variates using transformations with multiple roots, The
! American Statistician, 30(2), p. 88-90.

! This subroutine makes use of the subroutines:
!	- rnorm	: Box-Muller method for random normal draws

! Input arguments:
!	- mu		: mean parameter of the InvGaussian distribution
!	- lambda	: shape parameter of the InvGaussian distribution

! Output arguments:
!	- fn_val	: random InvGaussian variate

subroutine rinvgaus (mu, lambda, fn_val)

implicit none

! Precision statement
integer, parameter :: dp = kind(1.0d0)

! Input arguments:
real(dp), intent(in) :: mu, lambda

! Output arguments:
real(dp), intent(out) :: fn_val

! Internal arguments:
real(dp) :: nu, q, z


call rnorm(nu)

nu = nu*nu
q = mu + (nu*mu*mu)/(lambda*2.0_dp) - &
    mu/(2.0_dp*lambda)*sqrt(4.0_dp*mu*lambda*nu &
    + mu*mu*nu*nu)

call random_number(z)

if (z .le. (mu/(mu+q))) then
    fn_val = q
else
    fn_val = mu*mu/q
end if

end subroutine rinvgaus


!===========================================================================================


! Generates one random draw from the gamma distribution with
! mean = shape*scale. The algorithm is based on Marsaglia & Tsang 
! "A Simple Method for Gererating Gamma Variables" (2000)

! This subroutine makes use of the subroutines 
!	- rnor   	: generate one normal draw 

! Input arguments:
!	- shape		: shape parameter of the gamma distribution
!	- scale		: scale parameter of the gamma distribution

! Output arguments:
!	- fn_val	: random gamma variate Gamma(shape, scale)


subroutine rgamma (shape, scale, fn_val)

implicit none

! Precision statement:
integer, parameter :: dp = kind(1.0d0)

! Input arguments:
real(dp), intent(in) :: shape, scale

! Output arguments:
real(dp), intent(out) :: fn_val

! Internal arguments
real(dp) :: a, d, c, x, v, u
logical :: flag 


if (shape < 1.0_dp) then
  a = shape + 1.0_dp
else
  a = shape
end if

d = a - 1.0_dp/3.0_dp
c = 1.0_dp/sqrt(9.0_dp*d)

flag = .true.

do while (flag)
  v = 0.0_dp

  do while (v <= 0.0_dp)
    call rnorm(x)
    v = (1.0_dp + c*x)**3.0_dp
  end do

  call random_number(u)

  if (u < (1.0_dp-(0.0331_dp*(x**4.0_dp)))) then
    fn_val = d*v
    flag = .false.
  end if

  if (log(u) < ((0.5_dp*x*x) + (d*(1.0_dp - v + log(v))))) then
    fn_val = d*v
    flag = .false.
  end if

end do


if (shape < 1.0_dp) then
  call random_number(u)
  fn_val = (fn_val * (u**(1.0_dp/shape))) * scale
else
  fn_val = fn_val * scale
end if


end subroutine rgamma

end subroutine QRb_AL_mcmc
