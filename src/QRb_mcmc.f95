! Written by Dries F. Benoit 
! Faculty of economics and business administration
! Ghent University - BELGIUM
! LAST UPDATE (dd/mm/yy): 29/03/11

! MCMC sampler for binary quantile regression. This code is based on:
!! Benoit, D.F. and Van den Poel, D. (2011). Binary quantile
!! regression: A Bayesian approach based on the Asymmetric Laplace
!! distribution. Journal of Applied Econometrics (in press).

! Input arguments:
!	- n		: number of units of analysis
!	- nvar		: number of independent variables
!	- r		: number of MCMC iterations
!	- keep		: thinning parameter of MCMC draws
!	- y		: binary dependent variable
!	- p		: quantile of interest
!	- step		: Metropolis-Hastings stepsize for beta
!	- x		: matrix of regressors (1:n, 1:nvar)
!	- betabar	: prior mean for the regression parameters
!	- rooti		: solve(chol2inv(chol(prior precision matrix)))

! Output arguments:
!	- betadraw	: the matrix of regression parameter estimates
!	- loglike	: the vector of loglikelihoods over the subsequent draws
!	- rejrate	: the rejection rate of the MH-algorithm for beta


SUBROUTINE QRb_mcmc (n, nvar, r, keep, y, p, step, x, betabar, rooti, betadraw,loglike, rejrate)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)


! Input arguments:
INTEGER, INTENT(IN) :: n, r, nvar, keep
INTEGER, INTENT(IN), DIMENSION(1:n) :: y
REAL(dp), INTENT(IN) :: p, step
REAL(dp), INTENT(IN), DIMENSION(1:nvar) :: betabar
REAL(dp), INTENT(IN), DIMENSION(1:n,1:nvar) :: x
REAL(dp), INTENT(IN), DIMENSION(1:nvar,1:nvar) :: rooti

! Output arguments:
REAL(dp), INTENT(OUT) :: rejrate
REAL(dp), INTENT(OUT), DIMENSION(1:(r/keep)) :: loglike
REAL(dp), INTENT(OUT), DIMENSION(1:(r/keep),1:nvar) :: betadraw

! Internal arguments:
INTEGER :: i1, i2, naccept
REAL(dp) :: llold, u2, llnew, priornew, priorold, lldif, alpha, unif
REAL(dp), DIMENSION(1:n) :: z
REAL(dp), DIMENSION(1:nvar) :: oldbeta, u1, newbeta


! -- SET STARTING VALUES
! ---- beta
oldbeta = 0.0_dp

! ---- latent z
z = 0.0_dp

! ---- llold
CALL evalALD(n, z, 0.0_dp, 1.0_dp, p, llold)

! ---- naccept
naccept = 0


! -- START OF MCMC CHAIN
DO i1 = 1,r

  ! draw new value for beta
  DO i2 = 1,nvar
    CALL rnorm(u1(i2))
  END DO
  newbeta = oldbeta + step * u1
  
  CALL evalALD(n, (z - MATMUL(x, newbeta)), 0.0_dp, 1.0_dp, p, llnew)
  CALL logdmnorm(nvar, newbeta, betabar, rooti, priornew)
  CALL logdmnorm(nvar, oldbeta, betabar, rooti, priorold)

  lldif = llnew + priornew - llold - priorold

  ! accept or reject new draws
  alpha = MIN(1.0_dp, EXP(lldif))

  IF (alpha .LT. 1.0_dp) THEN
    CALL RANDOM_NUMBER(unif)
  ELSE
    unif = 0.0_dp
  END IF

  ! update variables if accepted
  IF (unif .LE. alpha) THEN
    llold = llnew
    oldbeta = newbeta
    naccept = naccept + 1
  END IF

  ! draw latent z
  DO i2 = 1,n
    CALL  rTALD (y(i2), DOT_PRODUCT(x(i2,1:nvar),oldbeta), 1.0_dp, p, z(i2))
  END DO

  CALL evalALD(n, (z - MATMUL(x, newbeta)), 0.0_dp, 1.0_dp, p, llold)

  IF (MOD(i1, keep) .EQ. 0) THEN
    betadraw((i1/keep),1:nvar) = oldbeta
    loglike(i1/keep) = llold
  END IF

END DO

rejrate = 1.0_dp - REAL(naccept,dp)/REAL(r,dp)


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

! Generate random number from the Truncated Asymmetric Laplace distribution
! Truncation depends on:
! 	IF (y==0) THEN (fn_val<0)
! 	IF (y==1) THEN (fn_val < 0)

! This subroutine makes use of the subroutines
!	pALD		: distribution function of the ALD
!	qALD		: quantile function of the ALD

! Input arguments:
!	- y		: defines truncation
!	- mu		: location parameter of the ALD
!	- sigma		: scale parameter of the ALD
!	- p		: quantile parameter of the ALD

! Output arguments:
!	- fn_val	: random number from truncated ALD



SUBROUTINE rTALD (y, mu, sigma, p, fn_val)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
INTEGER, INTENT(IN) :: y
REAL(dp), INTENT(IN) :: mu, sigma, p

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val

! Internal arguments:
REAL(dp) :: nulp, gr, min, max, u



nulp = -mu/sigma

IF ((y == 0) .AND. (nulp <= 0.0_dp)) THEN
  CALL RANDOM_NUMBER(u)
  fn_val = -(ABS(nulp) + (-LOG(u))/(1.0_dp-p))
ELSE IF ((y == 1) .AND. (nulp >= 0.0_dp)) THEN
  CALL RANDOM_NUMBER(u)
  fn_val = nulp + (-LOG(u))/p
ELSE
  CALL PALD(nulp, 0.0_dp, 1.0_dp, p, gr)
  IF (y == 1) THEN
    min = gr
    max = 1.0_dp
  ELSE
    min = 0.0_dp
    max = gr
  END IF

  CALL RANDOM_NUMBER(u)
  u = min + (max-min)*u
  CALL QALD(u, 0.0_dp, 1.0_dp, p, fn_val)
END IF

fn_val = mu + fn_val*sigma

END SUBROUTINE rTALD

!===========================================================================================

! This code evaluates the Asymmetric Laplace distribution function

! Input arguments:
!	- q		: location to be evaluated
!	- mu		: location parameter of the ALD
!	- sigma		: scale parameter of the ALD
!	- p		: quantile parameter of the ALD

! Output arguments:
!	- fn_val	: prob(-Inf, p] under ALD(mu, sigma, p)



SUBROUTINE pALD(q, mu, sigma, p, fn_val)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
REAL(dp), INTENT(IN) :: q, mu, sigma, p

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val



IF (q <= mu) THEN
  fn_val = p * EXP((1.0_dp-p) / sigma * (q - mu))
ELSE
  fn_val = 1.0_dp - ((1.0_dp-p) * EXP(-p / sigma * (q - mu)))
END IF

END SUBROUTINE pALD

!===========================================================================================

! This function evaluates the Asymmetric Laplace quantile function

! Input arguments:
!	- pp		: quantile to be evaluated
!	- mu		: location parameter of the ALD
!	- sigma		: scale parameter of the ALD
!	- p		: quantile parameter of the ALD

! Output arguments:
!	- fn_val	: evaluation of quantile function for ALD(mu, sigma, p)



SUBROUTINE qALD(pp, mu, sigma, p, fn_val)

IMPLICIT NONE

! Precision statement:
INTEGER, PARAMETER :: dp = KIND(1.0d0)

! Input arguments:
REAL(dp), INTENT(IN) :: pp, mu, sigma, p

! Output arguments:
REAL(dp), INTENT(OUT) :: fn_val


IF (pp <= p) THEN
  fn_val = mu + (sigma / (1.0_dp - p) * LOG(pp/p))
ELSE
  fn_val = mu - (sigma / p * LOG((1.0_dp - pp)/(1.0_dp - p)))
END IF

END SUBROUTINE qALD

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


END SUBROUTINE QRb_mcmc
