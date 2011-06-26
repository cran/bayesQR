! Written by 
! Dries F. Benoit 
! Faculty of economics and business administration
! Ghent University - BELGIUM
! and
! Rahim Al-Hamzawi
! Department of mathematics
! Brunel University - UK

! LAST UPDATE (dd/mm/yy): 25/06/11

! MCMC sampler for quantile regression with adaptive lasso. 
! This code is based on the following paper:
! Alhamzawi R, Yu K, Benoit DF. 2011. Bayesian adaptive LASSO quantile regression.
! Statistical Modeling (forthcoming).

! Input arguments:
!	- n		: number of units of analysis
!	- k		: number of independent variables
!	- r		: number of MCMC iterations
!	- keep		: thinning parameter of MCMC draws
!	- y		: dependent variable
!	- p		: quantile of interest
!	- step_delta	: Metropolis-Hastings stepsize for delta
!	- x		: matrix of regressors (1:n, 1:k)
!       - sig_shape     : shape hyperparameter for sigma
!       - sig_rate      : rate hyperparameter for sigma

! Output arguments:
!	- betadraw	: the matrix of regression parameter estimates
!       - lambda2draw   : the matrix of lambda2draws
!       - deltadraw     : the vector of delta's
!       - taudraw       : the vector of tau's
!       - rejrate       : the rejection rate of the M-H step

SUBROUTINE QRc_AL_mcmc (n, k, r, keep, y, p, step_delta, x, sig_shape, sig_rate, &
                        betadraw, lambda2draw, deltadraw, taudraw, rejrate)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
INTEGER, INTENT(IN) :: n, k, r, keep
REAL(dp), INTENT(IN) :: p, sig_shape, sig_rate, step_delta
REAL(dp), INTENT(IN), DIMENSION(1:n) :: y
REAL(dp), INTENT(IN), DIMENSION(1:n,1:k) :: x

! Output arguments:
REAL(dp), INTENT(OUT) :: rejrate
REAL(dp), INTENT(OUT), DIMENSION(1:(r/keep)) :: deltadraw, taudraw
REAL(dp), INTENT(OUT), DIMENSION(1:(r/keep),1:k) :: betadraw, lambda2draw

! Internal arguments:
INTEGER :: naccept, i1, i2, i3, cnt
INTEGER, DIMENSION(1:(k-1)) :: selec
REAL(dp) :: sigma, delta, tau, theta, phisq, quant1, tlambda1, tmu2, tlambda2, &
            tvar, tmean, rnor, tshape, trate1, deltanew, llnew, llold, lldiff, &
            alpha, unif
REAL(dp), DIMENSION(1:n) :: z, tmu1
REAL(dp), DIMENSION(1:k) :: beta, s, lambda2, trate2


! -- SET STARTING VALUES
z = 1.0_dp
beta = 1.0_dp
s = 1.0_dp
lambda2 = 1.0_dp
sigma = 1.0_dp
delta = 1.0_dp
tau = 1.0_dp
naccept = 0


! -- CALCULATE USEFUL QUANTITIES
theta = (1.0_dp - 2.0_dp*p)/(p*(1.0_dp-p))
phisq = 2.0_dp/(p*(1.0_dp-p))
quant1 = theta**2.0_dp + 2.0_dp*phisq

! -- START OF MCMC CHAIN
DO i1 = 1,r

  ! draw the latent z
  tmu1 = sqrt(quant1/(y-matmul(x, beta))**2.0_dp)
  tlambda1 = sigma*quant1/phisq
  DO i2 = 1,n
    CALL rInvGaussian(tmu1(i2), tlambda1, z(i2)) 
  END DO
  z = z**(-1.0_dp)

  ! draw mixing variable s
  DO i2 = 1,k
    tmu2 = sqrt(beta(i2)**2.0_dp*lambda2(i2)/sigma)
    IF (tmu2 .LT. 0.001_dp) THEN
      tmu2 = 0.001_dp
    END IF
    tlambda2 = beta(i2)**2.0_dp
    IF (tlambda2 .LT. 0.001_dp) THEN
      tlambda2 = 0.001_dp
    END IF
    CALL rInvGaussian(tmu2, tlambda2, s(i2))
  END DO

  ! draw the regression parameters beta
  DO i2 = 1,k
  cnt = 0
    DO i3 = 1,k
      IF (i2 .NE. i3) THEN
        cnt = cnt + 1
        selec(cnt) = i3
      END IF
    END DO

    tvar = (sigma*phisq**(-1.0_dp)*sum(x(1:n,i2)**2.0_dp*z**(-1.0_dp)) + &
           s(i2)**(-1.0_dp))**(-1.0_dp)
    tmean = tvar*sigma*phisq**(-1.0_dp)*sum(x(1:n,i2)*z**(-1.0_dp) * &
            (y - matmul(x(1:n,selec),beta(selec)) - theta*z))
    CALL rnorm(rnor)
    beta(i2) = tmean + rnor*sqrt(tvar)
  END DO

  ! draw sigma
  tshape = sig_shape + real(k,dp) + 3.0_dp/2.0_dp*real(n,dp)
  trate1 = sum(((y - matmul(x,beta) - theta*z)**2.0_dp/(2.0_dp*phisq*z)) + &
          z) + sum(s/(2.0_dp*lambda2)) + sig_rate
  CALL rgamma(tshape, trate1**(-1.0_dp), sigma)

  ! draw lambda2
  tshape = 1.0_dp + delta
  trate2 = s*sigma/2.0_dp + tau
  DO i2 = 1,k
    CALL rgamma(tshape, trate2(i2)**(-1.0_dp), lambda2(i2))
  END DO
  lambda2 = lambda2**(-1.0_dp)

  ! draw of hyperparameter tau
  tshape = delta*real(k,dp)
  trate1 = sum(lambda2)
  CALL rgamma(tshape, trate1**(-1.0_dp), tau)

  ! draw of shape hyperparameter delta: Metropolis-Hastings
  deltanew = -1.0_dp
  DO WHILE (deltanew .LE. 0.0_dp)
    CALL rnorm(rnor)
    deltanew = delta + rnor*step_delta
  END DO

  CALL LLdelta(deltanew, k, tau, lambda2, llnew)
  CALL LLdelta(delta, k, tau, lambda2, llold)

  lldiff = llnew - llold
  alpha = min(1.0_dp, exp(lldiff))

  IF (alpha .LT. 1.0_dp) THEN
    CALL random_number(unif)
  ELSE
    unif = 0.0_dp
  END IF

  IF (unif .LE. alpha) THEN
    delta = deltanew
    naccept = naccept + 1
  END IF

  ! write current draws to output arrays
  IF (mod(i1, keep) .EQ. 0) THEN
    betadraw((i1/keep),1:k) = beta
    lambda2draw((i1/keep),1:k) = lambda2
    deltadraw(i1/keep) = delta
    taudraw(i1/keep) = tau
  END IF

  rejrate = 1.0_dp - real(naccept,dp)/real(r,dp)

END DO

!=========================================================================

CONTAINS

!=========================================================================

! This code generates one draw from the standard normal 
! distribution. Note that more efficient code is possible
! when more than one normal draw is required.
! This code is based on the Box-Muller method.

! Output arguments:
!	- fn_val	: random draw from N(0,1) distribution

SUBROUTINE rnorm(fn_val)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Internal arguments:
REAL(dp) :: pi
REAL(dp), DIMENSION(1:2) :: u

pi = 3.14159265358979323846_dp

CALL random_number(u)

fn_val = sqrt(-2*log(u(1))) * cos(2*pi*u(2))

END SUBROUTINE rnorm

!=========================================================================

! This code generates one random draw from the inverse Gaussian distribution.
! The algorithm is based on: Michael, Schucany & Haas (1976), Generating
! random variates using transformations with multiple roots, The
! American Statistician, 30(2), p. 88-90.

! Input arguments:
!	- mu		: mean parameter of the InvGaussian distribution
!	- lambda	: shape parameter of the InvGaussian distribution

! Output arguments:
!	- fn_val	: random InvGaussian variate

SUBROUTINE rInvGaussian (mu, lambda, fn_val)

IMPLICIT NONE

! Precision statement
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
REAL(dp), INTENT(IN) :: mu, lambda

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Internal arguments:
REAL(dp) :: nu, q, z

CALL rnorm(nu)

nu = nu**2.0_dp
q = mu + (nu*mu**2.0_dp)/(lambda*2.0_dp) - &
    mu/(2.0_dp*lambda)*sqrt(4.0_dp*mu*lambda*nu + mu**2.0_dp*nu**2.0_dp)
nu = mu/(mu+q)

CALL random_number(z)

IF (z .LE. nu) THEN
    fn_val = q
ELSE
    fn_val = mu**2.0_dp/q
END IF

END SUBROUTINE rInvGaussian

!=========================================================================

! Implementation of the Lanczos approximation of the gamma
! function. Only use this code when the goal is to compute
! the LOGgamma, i.e. log(lancz_gamma(x)). Imprecise as 
! approximation for the gamma function for larger x.
! See:
! Lanczos, Cornelius (1964). "A Precision Approximation of the 
! Gamma Function". SIAM Journal on Numerical Analysis series B 
! (Society for Industrial and Applied Mathematics) 1: 86-96.

! Input arguments:
!   - x       : point to evaluate

! Output arguments:
!   - fn_val  : LOG of Lanczos approximation of the gamma function

RECURSIVE SUBROUTINE lancz_gamma(x, fn_val)

IMPLICIT NONE   

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
REAL(dp), INTENT(IN) :: x

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Local arguments:
INTEGER :: i1
REAL(dp) :: t, w, a, b
REAL(dp), DIMENSION(1:8) :: c
INTEGER, PARAMETER :: cg = 7
REAL(dp), PARAMETER :: pi = 3.14159265358979324_dp
REAL(dp), DIMENSION(0:8), PARAMETER :: p = &
        (/ 0.99999999999980993_dp, 676.5203681218851_dp, -1259.1392167224028_dp, &
        771.32342877765313_dp, -176.61502916214059_dp, 12.507343278686905_dp, &
        -0.13857109526572012_dp, 9.9843695780195716D-6, 1.5056327351493116D-7 /)

a = x

IF (a < .5_dp) THEN
        CALL lancz_gamma(1.0_dp - a, b)
        fn_val = pi / (sin(pi*a) * b)
ELSE
        a = a - 1.0_dp
        c(1) = a + 1.0_dp
        DO i1 = 1,7
          c(i1+1) = c(i1) + 1.0_dp
        END DO
        t = p(0) + sum(p(1:8)/c)
        w = a + REAL(cg,dp) + .5_dp
        fn_val = sqrt(2.0_dp*pi) * w**(a+.5_dp) * exp(-w) * t
END IF

END SUBROUTINE lancz_gamma

!=========================================================================

! This code evaluates the logliklihood of the delta parameter

SUBROUTINE LLdelta (x, k, tau, lambda2, fn_val)

IMPLICIT NONE   

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
INTEGER, INTENT(IN) :: k
REAL(dp), INTENT(IN) :: x, tau
REAL(dp), INTENT(IN), DIMENSION(1:k) :: lambda2

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Local arguments:
REAL(dp) :: w

CALL lancz_gamma(x, w)
fn_val = -real(k,dp)*log(w) + real(k,dp)*x*log(tau) - 2.0_dp*x* &
         sum(log(sqrt(lambda2)))

END SUBROUTINE LLdelta

!=========================================================================

! Generates one random draw from the gamma distribution with
! mean = shape*scale. The algorithm is based on Marsaglia & Tsang 
! "A Simple Method for Gererating Gamma Variables" (2000)

! Input arguments:
!	- shape		: shape parameter of the gamma distribution
!	- scale		: scale parameter of the gamma distribution

! Output arguments:
!	- fn_val	: random gamma variate Gamma(shape, scale)

SUBROUTINE rgamma (shape, scale, fn_val)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
REAL(dp), INTENT(IN) :: shape, scale

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Internal arguments
REAL(dp) :: a, d, c, x, v, u
LOGICAL :: flag 


IF (shape < 1.0_dp) THEN
  a = shape + 1.0_dp
ELSE
  a = shape
END IF

d = a - 1.0_dp/3.0_dp
c = 1.0_dp/SQRT(9.0_dp*d)

flag = .TRUE.

DO WHILE (flag)
  v = 0.0_dp

  DO WHILE (v <= 0.0_dp)
    CALL rnorm(x)
    v = (1.0_dp + c*x)**3.0_dp
  END DO

  CALL RANDOM_NUMBER(u)

  IF (u < (1.0_dp-(0.0331_dp*(x**4.0_dp)))) THEN
    fn_val = d*v
    flag = .FALSE.
  END IF

  IF (LOG(u) < ((0.5_dp*x*x) + (d*(1.0_dp - v + LOG(v))))) THEN
    fn_val = d*v
    flag = .FALSE.
  END IF

END DO


IF (shape < 1.0_dp) THEN
  CALL RANDOM_NUMBER(u)
  fn_val = (fn_val * (u**(1.0_dp/shape))) * scale
ELSE
  fn_val = fn_val * scale
END IF


END SUBROUTINE rgamma

!=========================================================================

END SUBROUTINE QRc_AL_mcmc

