! Written by Dries F. Benoit 
! Faculty of economics and business administration
! Ghent University - BELGIUM
! LAST UPDATE (dd/mm/yy): 29/03/11

! MCMC sampler for quantile regression. This code is based on:
! Yu K, Moyeed RA. 2001. Bayesian Quantile Regression. Statistics 
! and Probability Letters 54(4): 437-447.

! Input arguments:
!	- n		: number of units of analysis
!	- nvar		: number of independent variables
!	- r		: number of MCMC iterations
!	- keep		: thinning parameter of MCMC draws
!	- y		: dependent variable
!	- p		: quantile of interest
!	- step1		: Metropolis-Hastings stepsize for beta
!	- step2		: Metropolis-Hastings stepsize for sigma
!	- x		: matrix of regressors (1:n, 1:nvar)
!	- betabar	: prior mean for the regression parameters
!	- rooti		: solve(chol2inv(chol(precision matrix)))
!	- nu		: prior invChiSq d.f. parameter
!	- ssq		: prior invChiSq scale parameter

! Output arguments:
!	- betadraw	: the matrix of regression parameter estimates
!	- sigdraw	: the vector of scale estimates
!	- loglike	: the vector of loglikelihoods over the subsequent draws
!	- rejrate1	: the rejection rate of the MH-algorithm for beta
!	- rejrate2	: the rejection rate of the MH-algorithm for sigma


SUBROUTINE QRc_mcmc (n, nvar, r, keep, y, p, step1, step2, x, betabar, rooti, &
                     nu, ssq, betadraw, sigdraw, loglike, rejrate1, rejrate2)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)


! Input arguments:
INTEGER, INTENT(IN) :: n, r, nvar, keep
REAL(dp), INTENT(IN) :: p, step1, step2, nu, ssq
REAL(dp), INTENT(IN), DIMENSION(1:n) :: y
REAL(dp), INTENT(IN), DIMENSION(1:nvar) :: betabar
REAL(dp), INTENT(IN), DIMENSION(1:n,1:nvar) :: x
REAL(dp), INTENT(IN), DIMENSION(1:nvar,1:nvar) :: rooti

! Output arguments:
REAL(dp), INTENT(OUT) :: rejrate1, rejrate2
REAL(dp), INTENT(OUT), DIMENSION(1:(r/keep)) :: sigdraw, loglike
REAL(dp), INTENT(OUT), DIMENSION(1:(r/keep),1:nvar) :: betadraw

! Internal arguments:
INTEGER :: i1, i2, naccept1, naccept2
REAL(dp) :: oldsig, llold, u2, newsig, llnew, priornew, priorold, lldif, alpha, unif
REAL(dp), DIMENSION(1:nvar) :: oldbeta, u1, newbeta


! -- SET STARTING VALUES
! ---- beta
oldbeta = 0.0_dp

! ---- sigma
oldsig = 1.0_dp

! ---- llold
CALL evalALD(n, y, 0.0_dp, oldsig, p, llold)

! ---- naccept
naccept1 = 0
naccept2 = 0


! -- START OF MCMC CHAIN
DO i1 = 1,r

  ! draw new value for beta
  DO i2 = 1,nvar
    CALL rnorm(u1(i2))
  END DO
  newbeta = oldbeta + step1 * u1
  
  CALL evalALD(n, (y - MATMUL(x, newbeta)), 0.0_dp, oldsig, p, llnew)
  CALL logdmnorm(nvar, newbeta, betabar, rooti, priornew)
  CALL logdmnorm(nvar, oldbeta, betabar, rooti, priorold)

  lldif = llnew + priornew - llold - priorold

  ! accept or reject new draws
  alpha = min(1.0_dp, exp(lldif))

  IF (alpha .LT. 1.0_dp) THEN
    CALL RANDOM_NUMBER(unif)
  ELSE
    unif = 0.0_dp
  END IF

  ! update variables if accepted
  IF (unif .LE. alpha) THEN
    llold = llnew
    oldbeta = newbeta
    naccept1 = naccept1 + 1
  END IF


  ! draw new value for sigma
  newsig = -1.0_dp
  DO WHILE (newsig .LE. 0.0_dp) 
    CALL rnorm(u2)
    newsig = oldsig + step2 * u2
  END DO

  CALL evalALD(n, (y - MATMUL(x, oldbeta)), 0.0_dp, newsig, p, llnew)
  CALL logdichisq(newsig, nu, ssq, priornew)
  CALL logdichisq(oldsig, nu, ssq, priorold)

  lldif = llnew + priornew - llold - priorold

  ! accept or reject new draws
  alpha = min(1.0_dp, exp(lldif))

  IF (alpha .LT. 1.0_dp) THEN
    CALL RANDOM_NUMBER(unif)
  ELSE
    unif = 0.0_dp
  END IF

  ! update variables if accepted
  IF (unif .LE. alpha) THEN
    llold = llnew
    oldsig = newsig
    naccept2 = naccept2 + 1
  END IF

  IF (MOD(i1, keep) .EQ. 0) THEN
    betadraw((i1/keep),1:nvar) = oldbeta
    sigdraw(i1/keep) = oldsig
    loglike(i1/keep) = llold
  END IF

END DO

rejrate1 = 1.0_dp - REAL(naccept1,dp)/REAL(r,dp)
rejrate2 = 1.0_dp - REAL(naccept2,dp)/REAL(r,dp)


!===========================================================================================

contains

!===========================================================================================


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

!===========================================================================================

! Evaluate the likelihood of the ALD
! Function returns the log likelihood

! Input arguments:
!	- n		: number of points to evaluate
!	- u		: points to be evaluated
!	- mu		: location parameter of the ALD
!	- sigma		: scale parameter of the ALD
!	- p		: quantile parameter of the ALD

! Output arguments:
!	- fn_val	: LL of the ALD evaluated at the u's

! Note: expected input sizes are programmed oddly!

SUBROUTINE evalALD(n, u, mu, sigma, p, fn_val)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
INTEGER, INTENT(IN) :: n
REAL(dp), INTENT(IN) :: mu, sigma, p
REAL(dp), INTENT(IN), DIMENSION(1:n) :: u

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Internal arguments:
REAL(dp), DIMENSION(1:n) :: u1, rho



u1 = (u-mu)/sigma
rho = (ABS(u1) + (2.0_dp*p-1.0_dp)*u1)/2.0_dp

fn_val = SUM(LOG(((p*(1.0_dp-p))/sigma) *EXP(-rho)))

END SUBROUTINE evalALD

!===========================================================================================

! Computes the log of a multi-variate normal density
! Note: this is not vectorized input!!

! Input arguments:
! 	- k	:	number of dimensions
! 	- x	:	point to evaluate
! 	- mu	:	mean vector of mvn
! 	- rooti	:	solve(chol2inv(chol(precision matrix)))

! Output arguments:
!	- fn_val:	log of MVN(mu, rooti) evaluated at x


SUBROUTINE logdmnorm (k, x, mu, rooti, fn_val)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
INTEGER, INTENT(IN) :: k
REAL(dp), INTENT(IN), DIMENSION(1:k) :: x, mu
REAL(dp), INTENT(IN), DIMENSION(1:k,1:k) :: rooti

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Internal arguments:
INTEGER :: i1
REAL(dp) :: pi
REAL(dp), DIMENSION(1:k) :: z, diag



pi = 3.14159265358979323846_dp

! Define the diagonal of rooti
DO i1 = 1,k
  diag(i1) = rooti(i1,i1)
END DO

! Compute log of MVN
z = MATMUL(TRANSPOSE(rooti), (x-mu))
fn_val = REAL(-k,dp)/2.0_dp * LOG(2.0_dp*pi) - 0.5_dp * DOT_PRODUCT(z,z) + SUM(LOG(diag))

END SUBROUTINE logdmnorm

!===========================================================================================

! Computes the log of a inverted Chi-square distribution

! Input arguments:
! 	- x	:	point to evaluate
! 	- nu	:	degrees of freedom parameter
! 	- ssq	:	scale parameter

! Output arguments:
!	- fn_val:	log of invChisq(nu, ssq) evaluated at x


SUBROUTINE logdichisq (x, nu, ssq, fn_val)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
REAL(dp), INTENT(IN) :: x, nu, ssq

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Internal arguments:
REAL(dp) :: a, b

a = -LOG_GAMMA(nu/2.0_dp) + (nu/2.0_dp) * LOG((nu * ssq)/2.0_dp)
b = -((nu/2.0_dp) + 1.0_dp) * LOG(x) - (nu * ssq)/(2.0_dp * x)
fn_val = a + b


END SUBROUTINE logdichisq

!===========================================================================================


END SUBROUTINE QRc_mcmc
