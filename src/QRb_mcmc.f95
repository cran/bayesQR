! Written by Dries F. Benoit 
! Faculty of economics and business administration
! Ghent University - BELGIUM

! MCMC sampler for quantile regression. This code is and improved
! version (Gibbs instead of Metropolis-Hastings) of the mcmc sampler
! proposed in:
! Yu K, Moyeed RA. 2001. Bayesian Quantile Regression. Statistics 
! and Probability Letters 54(4): 437-447.

! Input arguments:
!	- n		      : number of units of analysis
!	- k   		  : number of independent variables
!	- r		      : number of MCMC iterations
!	- keep		  : thinning parameter of MCMC draws
!	- y		      : dependent variable
!	- p		      : quantile of interest
!	- x		      : matrix of regressors (1:n, 1:nvar)
!	- beta0    	: prior mean for the regression parameters
!	- V0i 	    : prior inverse covariance matrix for regression parameters

! Output arguments:
!	- betadraw	: the matrix of regression parameter estimates


subroutine QRb_mcmc (n, k, r, keep, y, p, x, beta0, V0i, betadraw)

implicit none

! Precision statement:
integer, parameter :: dp = kind(1.0d0)

! Input arguments:
integer, intent(in) :: n, r, k, keep
integer, intent(in), dimension(n) :: y
real(dp), intent(in) :: p
real(dp), intent(in), dimension(k) :: beta0
real(dp), intent(in), dimension(n,k) :: x
real(dp), intent(in), dimension(k,k) :: V0i

! Output arguments:
real(dp), intent(out), dimension(r/keep,k) :: betadraw

! Internal arguments:
integer :: i1, i2, i3, ok
real(dp) :: theta, tausq, lambda
real(dp), dimension(k) :: beta, betabar
real(dp), dimension(n) :: v, ystar, Xbeta, mu
real(dp), dimension(k,k) :: Vbaro1 


! -- USEFUL QUANTITIES 
theta = (1.0_dp - 2.0_dp*p)/(p*(1.0_dp - p))
tausq = 2.0_dp/(p*(1.0_dp - p))

! -- SET STARTING VALUES
beta = 0.0_dp
v = 1.0_dp
ystar = 0.0_dp

! -- START OF MCMC CHAIN
do i1 = 1,r
  ! Simulate new values for ystar 
  !ooooooooooooooooooooooooooooooo
    Xbeta = matmul(X,beta)
    do i2 = 1,n
      if (y(i2) == 0) then
        call rtnorm_geweke(0.0_dp, .false., Xbeta(i2)+theta*v(i2), sqrt(tausq)&
                           *sqrt(v(i2)), ystar(i2))
      else
        call rtnorm_geweke(0.0_dp, .true., Xbeta(i2)+theta*v(i2), sqrt(tausq)&
                           *sqrt(v(i2)), ystar(i2))
      endif
    enddo
    

  ! Simulate new values for v 
  !ooooooooooooooooooooooooooo
    lambda = 2.0_dp + (theta**2.0_dp)/tausq
    mu = sqrt(lambda/((ystar - matmul(X,beta))**2.0_dp/tausq))
    do i2 = 1,n
      call rinvgaus(mu(i2),lambda,v(i2))
    enddo
    v = 1.0_dp/v

  ! Simulate new value for beta
  !ooooooooooooooooooooooooooooo
    Vbaro1 = 0.0_dp
    do i2 = 1,n
      ! Code below: outer/tensor product (see Numerical Recipes in Fortran 90)
      Vbaro1 = Vbaro1 + spread(X(i2,:),dim=2,ncopies=k)*spread(X(i2,:),dim=1,&
                        ncopies=k)/(tausq*v(i2))
    enddo
    Vbaro1 = Vbaro1 + V0i
  
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
      betabar = betabar + X(i2,:)*(ystar(i2)-theta*v(i2))/(tausq*v(i2))
    enddo
    betabar = betabar + matmul(V0i,beta0)
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

end subroutine QRb_mcmc
