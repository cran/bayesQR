QRb.AL <- function(Data, Prior, Mcmc){

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
		cat("length y: ",n, "ncol(X): ", ncol(X))
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
	c = 0.01
	d = 0.01
    }
    else {
        if (is.null(Prior$lambdasq_shape)) {
	    c = 0.01
        }
        else {
            c = Prior$lambdasq_shape
        }
        if (is.null(Prior$lambdasq_scale)) {
	    d = 1/0.01
        }
        else {
            d = 1/Prior$lambdasq_scale
        }
    }

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
    x <- as.double(X)
    c <- as.double(c)
    d <- as.double(d)
    betadraw <- double(nvar*r/keep)

    ## Call Fortran routine
    fn_val <- .Fortran("QRb_AL_mcmc", n, nvar, r, keep, y, p, x, c, d, betadraw) 

    return(list(method="QRb.AL",
		            p=p,
		            betadraw=matrix(fn_val[[10]], nrow=r/keep, ncol=nvar)
		            )
					)
}
