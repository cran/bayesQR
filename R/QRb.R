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
    y = as.vector(Data$y)
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
        beta0 = rep(0, nvar)
        V0 = 100 * diag(nvar)
    }
    else {
        if (is.null(Prior$beta0)) {
            beta0 = rep(0, nvar)
        }
        else {
            beta0 = Prior$beta0
        }
        if (is.null(Prior$V0)) {
            V0 = 100 * diag(nvar)
        }
        else {
            V0 = Prior$V0
        }
    }
    if (ncol(V0) != nrow(V0) || ncol(V0) != nvar || nrow(V0) != nvar) {
        pandterm(paste("bad dimensions for V0", dim(V0)))
    }
    if (length(beta0) != nvar) {
        pandterm(paste("beta0 wrong length, length= ", length(beta0)))
    }
    V0i = chol2inv(chol(V0))

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
    }


# Start of algorithm:

    ## Assign correct variable types
    n <- as.integer(n)
    nvar <- as.integer(nvar)
    r <- as.integer(r)
    keep <- as.integer(keep)
    y <- as.integer(y)
    p <- as.double(p)
    X <- as.double(X)
    beta0 <- as.double(beta0)
    V0i <- as.double(V0i)
    betadraw <- double(nvar*r/keep)

    ## Call Fortran routine
    fn_val <- .Fortran("QRb_mcmc", n, nvar, r, keep, y, p, X, beta0, V0i, betadraw)


    return(list(method="QRb",
		            p=p,
		            betadraw=matrix(fn_val[[10]], nrow=r/keep, ncol=nvar)
		            )
					)
}
