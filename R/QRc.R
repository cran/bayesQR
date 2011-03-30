QRc <- function(Data, Prior, Mcmc){

# Error handling:

    pandterm = function(message) {
        stop(message, call. = FALSE)
    }
    if (missing(Data)) {
        pandterm("Requires Data argument -- list of y, X and p")
    }
    if (is.null(Data$X)) {
        pandterm("Requires Data element X")
    }
    X = Data$X
    if (is.null(Data$y)) {
        pandterm("Requires Data element y")
    }
    y = Data$y
    nvar = ncol(X)
    n = length(y)
    if (n != nrow(X)) {
        pandterm("length(y) ne nrow(X)")
    }
    if (is.null(Data$p)) {
        pandterm("Requires Data element p")
    }
    p = Data$p
    if ((p<=0)|(p>=1)) {
        pandterm("Unvalid quantile p")
    }

    if (missing(Prior)) {
        betabar = c(rep(0, nvar))
        A = 0.01 * diag(nvar)
        nu = 3
        ssq = var(y)
    }
    else {
        if (is.null(Prior$betabar)) {
            betabar = c(rep(0, nvar))
        }
        else {
            betabar = Prior$betabar
        }
        if (is.null(Prior$A)) {
            A = 0.01 * diag(nvar)
        }
        else {
            A = Prior$A
        }
        if (is.null(Prior$nu)) {
            nu = 3
        }
        else {
            nu = Prior$nu
        }
        if (is.null(Prior$ssq)) {
            ssq = var(y)
        }
        else {
            ssq = Prior$ssq
        }
    }
    if (ncol(A) != nrow(A) || ncol(A) != nvar || nrow(A) != nvar) {
        pandterm(paste("bad dimensions for A", dim(A)))
    }
    if (length(betabar) != nvar) {
        pandterm(paste("betabar wrong length, length= ", length(betabar)))
    }
    Ai = chol2inv(chol(A))
    rooti = solve(chol(Ai))

    if (missing(Mcmc)) {
        pandterm("Requires Mcmc argument")
    }
    else {
        if (is.null(Mcmc$R)) {
            pandterm("Requires Mcmc element R")
        }
        else {
            r = Mcmc$R
        }
        if (is.null(Mcmc$keep)) {
            keep = 1
        }
        else {
            keep = Mcmc$keep
        }
        if (is.null(Mcmc$step_beta)) {
            pandterm("Requires Mcmc element step_beta")
        }
        else {
            step1 = Mcmc$step_beta
        }
        if (is.null(Mcmc$step_sigma)) {
            pandterm("Requires Mcmc element step_sigma")
        }
        else {
            step2 = Mcmc$step_sigma
        }
    }


# Start of algorithm:

    ## Assign correct variable types
    n <- as.integer(n)
    nvar <- as.integer(nvar)
    r <- as.integer(r)
    keep <- as.integer(keep)
    y <- as.double(y)
    p <- as.double(p)
    step1 <- as.double(step1)
    step2 <- as.double(step2)
    X <- as.double(X)
    betabar <- as.double(betabar)
    rooti <- as.double(rooti)
    nu <- as.double(nu)
    ssq <- as.double(ssq)
    betadraw <- double(nvar*r/keep)
    sigdraw <- double(r/keep)
    loglike <- double(r/keep)
    rejrate1 <- double(1)
    rejrate2 <- double(1)

    ## Call FORTRAN routine
    fn_val <- .Fortran("QRc_mcmc", n, nvar, r, keep, y, p, step1, step2, X, betabar, 
			rooti, nu, ssq, betadraw, sigdraw, loglike, rejrate1, rejrate2)


    return(list(betadraw=matrix(fn_val[[14]], nrow=r/keep, ncol=nvar), 
			sigmadraw=fn_val[[15]], loglike=fn_val[[16]], 
			rejrate_beta=fn_val[[17]],rejrate_sigma=fn_val[[18]]))
}
