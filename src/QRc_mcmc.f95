! Written by Dries F. Benoit 
! Faculty of economics and business administration
! Ghent University - BELGIUM

! MCMC sampler for quantile regression. This code is and improved
! version (Gibbs instead of Metropolis-Hastings) of the mcmc sampler
! proposed in:
! Yu K, Moyeed RA. 2001. Bayesian Quantile Regression. Statistics 
! and Probability Letters 54(4): 437-447.

! Input arguments:
!	- n	    : number of units of analysis
!	- k   	    : number of independent variables
!	- r	    : number of MCMC iterations
!	- keep	    : thinning parameter of MCMC draws
!	- y	    : dependent variable
!	- p	    : quantile of interest
!	- x	    : matrix of regressors (1:n, 1:nvar)
!	- beta0     : prior mean for the regression parameters
!	- V0 	    : prior covariance matrix for regression parameters
!       - shape0    : prior invGamma shape parameter for sigma
!       - scale0    : prior invGamma scale parameter for sigma
!       - normal    : indicator if normal approximation is requested

! Output arguments:
!	- betadraw	: the matrix of regression parameter estimates
!	- sigmadraw	: the vector of scale estimates


subroutine QRc_mcmc (n, k, r, keep, y, p, x, beta0, V01, shape0, scale0, normal, betadraw, sigmadraw)

implicit none

! Precision statement:
integer, parameter :: dp = kind(1.0d0)

! Input arguments:
logical, intent(in) :: normal 
integer, intent(in) :: n, r, k, keep
real(dp), intent(in) :: p, shape0, scale0 
real(dp), intent(in), dimension(n) :: y
real(dp), intent(in), dimension(k) :: beta0
real(dp), intent(in), dimension(n,k) :: x
real(dp), intent(in), dimension(k,k) :: V01

! Output arguments:
real(dp), intent(out), dimension(r/keep) :: sigmadraw
real(dp), intent(out), dimension(r/keep,k) :: betadraw

! Internal arguments:
integer :: i1, i2, i3, ok
real(dp) :: theta, tausq, sigma, lambda, shapebar, scalebar
real(dp), dimension(k) :: beta, betabar
real(dp), dimension(n) :: v, mu 
real(dp), dimension(k,k) :: Vbaro1 


! -- USEFUL QUANTITIES 
theta = (1.0_dp - 2.0_dp*p)/(p*(1.0_dp - p))
tausq = 2.0_dp/(p*(1.0_dp - p))

! -- SET STARTING VALUES
beta = 0.0_dp
v = 1.0_dp
sigma = 1.0_dp

! -- START OF MCMC CHAIN
do i1 = 1,r

  ! Simulate new values for v 
  !ooooooooooooooooooooooooooo
    lambda = 2.0_dp/sigma + (theta**2.0_dp)/(tausq*sigma)
    mu = sqrt(lambda/((y - matmul(X,beta))**2.0_dp/(tausq*sigma)))
    do i2 = 1,n
      call rinvgaus(mu(i2),lambda,v(i2))
    enddo
    v = 1.0_dp/v

  ! Simulate new value for beta
  !ooooooooooooooooooooooooooooo
    Vbaro1 = 0.0_dp
    do i2 = 1,n
      ! Code below: outer/tensor product (see Numerical Recipes in Fortran 90)
      Vbaro1 = Vbaro1 + spread(X(i2,:),dim=2,ncopies=k)*spread(X(i2,:),dim=1,ncopies=k)/&
                        (tausq*sigma*v(i2))
    enddo
    Vbaro1 = Vbaro1 + V01
  
    ! Invert Vbaro1
    call dpotrf('U',k,Vbaro1,k,ok)
    call dpotri('U',k,Vbaro1,k,ok)
    do i2 = 1,(k-1)
      do i3 = (i2+1),k
        Vbaro1(i3,i2) = Vbaro1(i2,i3)
      enddo
    enddo
  
    betabar = 0.0_dp
    do i2 = 1,n
      betabar = betabar + X(i2,:)*(y(i2)-theta*v(i2))/(tausq*sigma*v(i2))
    enddo
    betabar = betabar + matmul(V01,beta0)
    betabar = matmul(Vbaro1,betabar)

    ! Cholesky factorization of Vbaro1
    call dpotrf('U',k,Vbaro1,k,ok)
    do i2 = 1,(k-1)
      do i3 = (i2+1),k
        Vbaro1(i3,i2) = 0.0_dp 
      enddo
    enddo

    do i2 = 1,k
      call rnorm(beta(i2))
    enddo

    beta = betabar + matmul(beta,Vbaro1)

  ! Simulate new value for sigma
  !oooooooooooooooooooooooooooooo
  if (.not. normal) then
    shapebar = (shape0*2.0_dp + 3.0_dp*real(n,dp))/2.0_dp
    scalebar = (scale0*2.0_dp + sum(v)*2.0_dp + sum((y - matmul(X,beta) -&
                theta*v)**2.0_dp/(tausq*v)))/2.0_dp
    if (scalebar < 1.0d-2) scalebar = 1.0d-2 !numerical stability
    if (scalebar > 1.0d2) scalebar = 1.0d2 !numerical stability
    call rgamma(shapebar,1.0_dp/scalebar,sigma)
    sigma = 1.0_dp/sigma
  endif

  ! Save current draw 
  !ooooooooooooooooooo
  if (mod(i1, keep) == 0) then
    betadraw((i1/keep),1:k) = beta
    sigmadraw(i1/keep) = sigma
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

end subroutine QRc_mcmc
