QRc.AL <- function(Data, Prior, Mcmc) {

# Error handling:

    pandterm = function(message) {
        stop(message, call. = FALSE)
    }
    if (missing(Data)) {
        pandterm("Requires Data argument -- list of y, X and p")
    } else {
        if (is.null(Data$X)) {
            pandterm("Requires Data element X")
        } else {
            x = Data$X
            k = ncol(x)
        }
        if (is.null(Data$y)) {
            pandterm("Requires Data element y")
        } else {
            y = Data$y
            n = length(y)
        }
        if (n != nrow(x)) {
            pandterm("length(y) ne nrow(X)")
        }
        if (is.null(Data$p)) {
            pandterm("Requires Data element p")
        } else {
            p = Data$p
        }
        if ((p<=0)|(p>=1)) {
            pandterm("Unvalid quantile p")
	}
    }

    if (missing(Prior)) {
        sig_shape = .01
        sig_rate = .01
    } else {
        if (is.null(Prior$sig_shape)) {
            sig_shape = .01
        } else {
            sig_shape = Prior$sig_shape
        }
        if (is.null(Prior$sig_rate)) {
            sig_rate = .01
        } else {
            sig_rate = Prior$sig_rate
        }
    }

    if (missing(Mcmc)) {
        pandterm("Requires Mcmc argument")
    } else {
        if (is.null(Mcmc$R)) {
            pandterm("Requires Mcmc element R")
        } else {
            r = Mcmc$R
        }
        if (is.null(Mcmc$keep)) {
            keep = 1
        } else {
            keep = Mcmc$keep
        }
        if (is.null(Mcmc$step_delta)) {
            pandterm("Requires Mcmc element step_beta")
        } else {
            step_delta = Mcmc$step_delta
        }
    }

# Start of algorithm

    ## Assign correct variable types
    n <- as.integer(n)
    k <- as.integer(k)
    r <- as.integer(r)
    keep <- as.integer(keep)
    y <- as.double(y)
    p <- as.double(p)
    step_delta <- as.double(step_delta)
    x <- as.double(x)
    sig_shape <- as.double(sig_shape)
    sig_rate <- as.double(sig_rate)
    betadraw <- double(k*r/keep)
    lambda2draw <- double(k*r/keep)
    deltadraw <- double(r/keep)
    taudraw <- double(r/keep)
    rejrate <- double(1)


    ## Call Fortran routine
    fn_val <- .Fortran("QRc_AL_mcmc", n, k, r, keep, y, p, step_delta, x, 
                        sig_shape, sig_rate, betadraw, lambda2draw, 
                        deltadraw, taudraw, rejrate)

    return(list(betadraw=matrix(fn_val[[11]], nrow=r/keep, ncol=k),
                lambda2draw=matrix(fn_val[[12]], nrow=r/keep, ncol=k),
                deltadraw=fn_val[[13]], taudraw=fn_val[[14]],
                rejrate=fn_val[[15]]))

}
