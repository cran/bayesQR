QRb <- function(Data, Prior, Mcmc){

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
    if ((sort(unique(y))[1]!=0)|(sort(unique(y))[2]!=1)) {
        pandterm("Unvalid dependent variable y")
    }
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
        if (is.null(Mcmc$step)) {
            pandterm("Requires Mcmc element step")
        }
        else {
            step = Mcmc$step
        }
    }


# Start of algorithm:

    ## Assign correct variable types
    n <- as.integer(n)
    nvar <- as.integer(nvar)
    r <- as.integer(r)
    keep <- as.integer(keep)
    y <- as.integer(y)
    p <- as.double(p)
    step <- as.double(step)
    X <- as.double(X)
    betabar <- as.double(betabar)
    rooti <- as.double(rooti)
    betadraw <- double(nvar*r/keep)
    loglike <- double(r/keep)
    rejrate <- double(1)

    ## Call FORTRAN routine
    fn_val <- .Fortran("QRb_mcmc", n, nvar, r, keep, y, p, step, X, betabar, 
                                   rooti, betadraw, loglike, rejrate)


    return(list(betadraw=matrix(fn_val[[11]], nrow=r/keep, ncol=nvar), 
                loglike=fn_val[[12]], rejrate=fn_val[[13]]))
}
