! Written by Dries F. Benoit 
! Faculty of economics and business administration
! Ghent University - BELGIUM

! MCMC sampler for quantile regression with lasso variable selection.
! This code is and  

! Input arguments:
!    - n          : number of units of analysis
!    - k          : number of independent variables
!    - r          : number of MCMC iterations
!    - keep       : thinning parameter of MCMC draws
!    - y          : dependent variable
!    - p          : quantile of interest
!    - x          : matrix of regressors (1:n, 1:nvar)
!    - a          : shape parameter of inverse gamma prior on sigma
!    - b          : scale parameter of inverse gamma prior on sigma
!    - c          : shape parameter of gamma prior on etasquared 
!    - d          : rate parameter of gamma prior on etasquared 
!    - normal     : indicator if normal approximation is requested 

! Output arguments:
!    - betadraw	  : the matrix of regression parameter estimates
!    - sigmadraw  : the vector of scale estimates


subroutine QRc_AL_mcmc (n, k, r, keep, y, p, x, a, b, c, d, normal, betadraw, sigmadraw)

implicit none

! Precision statement:
integer, parameter :: dp = kind(1.0d0)

! Input arguments:
logical, intent(in) :: normal
integer, intent(in) :: n, k, r, keep
real(dp), intent(in) :: p, a, b, c, d 
real(dp), intent(in), dimension(n) :: y
real(dp), intent(in), dimension(n,k) :: x

! Output arguments:
real(dp), intent(out), dimension(r/keep) :: sigmadraw 
real(dp), intent(out), dimension(r/keep,k) :: betadraw

! Internal arguments:
integer :: i1, i2
real(dp) :: theta, omegasq, sigma, sigsq, shape, scale
real(dp), dimension(k) :: s, beta, etasq, lambda, beta0
real(dp), dimension(n) :: nu, mu, ytilde 


! -- USEFUL QUANTITIES 
theta = (1.0_dp - 2.0_dp*p)/(p*(1.0_dp - p))
omegasq = 2.0_dp/(p*(1.0_dp - p))

! -- SET STARTING VALUES
nu = 1.0_dp
s = 1.0_dp
beta = 1.0_dp
sigma = 1.0_dp
etasq = 1.0_dp

! -- START OF MCMC CHAIN
do i1 = 1,r

  ! Simulate new values for nu 
  !oooooooooooooooooooooooooooo
    lambda(1) = theta**2/(sigma*omegasq) + 2.0_dp/sigma
    mu = sqrt(lambda(1)/((y-matmul(x,beta))**2/(sigma*omegasq)))
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
        ytilde = y - theta*nu - matmul(x,beta0)
        sigsq = 1.0_dp/(sum((x(:,i2)**2)/nu)/(sigma*omegasq) + 1.0_dp/s(i2))
        mu(1) = sigsq/(sigma*omegasq)*sum(ytilde*x(:,i2)/nu)
        call rnorm(beta(i2))
        beta(i2) = mu(1) + beta(i2)*sqrt(sigsq)
    enddo

  ! Simulate new values for sigma 
  !ooooooooooooooooooooooooooooooo
  if (.not. normal) then
    shape = a + 1.5_dp*real(n,dp)
    scale = sum(((y-matmul(X,beta)-theta*nu)**2)/(2.0_dp*omegasq*nu) + nu) + b
    if (scale < 1.0d-2) scale = 1.0d-2 !numerical stability
    if (scale > 1.0d2) scale = 1.0d2 !numerical stability
    call rgamma(shape,1.0_dp/scale,sigma)
    sigma = 1.0_dp/sigma
  endif

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
    sigmadraw((i1/keep)) = sigma 
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


! This code generates one draw from the standard uniform 
! distribution. It calls R's internal random number
! generation routine.

! Output arguments:
!	- fn_val	: random draw from U(0,1) distribution

subroutine runif(fn_val)

implicit none

! Precision statement:
integer, parameter :: dp = kind(1.0d0)

! Output arguments:
real(dp), intent(out) :: fn_val

! Internal arguments:
real(dp) :: unifrnd

call rndstart()
fn_val = unifrnd()
call rndend()

end subroutine runif



!===========================================================================================


! This code generates one draw from the standard Gaussian 
! distribution. It calls R's internal random number
! generation routine.

! Output arguments:
!	- fn_val	: random draw from N(0,1) distribution

subroutine rnorm(fn_val)

implicit none

! Precision statement:
integer, parameter :: dp = kind(1.0d0)

! Output arguments:
real(dp), intent(out) :: fn_val

! Internal arguments:
real(dp) :: normrnd

call rndstart()
fn_val = normrnd()
call rndend()

end subroutine rnorm


!===========================================================================================


! This code generates one random draw from the inverse Gaussian distribution.
! The algorithm is based on: Michael, Schucany & Haas (1976), Generating
! random variates using transformations with multiple roots, The
! American Statistician, 30(2), p. 88-90.

! This subroutine makes use of the subroutines:
!	- rnorm	: R's internal random normal draw generator

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

call runif(z)

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

  call runif(u)

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
  call runif(u)
  fn_val = (fn_val * (u**(1.0_dp/shape))) * scale
else
  fn_val = fn_val * scale
end if

end subroutine rgamma

end subroutine QRc_AL_mcmc
