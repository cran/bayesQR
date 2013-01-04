QRsummary <- function(QRobj,burnin=1,credint=c(0.025,0.975),traceplot=FALSE){

    # Error handling
    pandterm <- function(message) {
        stop(message, call. = FALSE)
    }
    if (missing(QRobj)){
        pandterm("Requires QRobj as provided by 'QRc', 'QRc.AL', 'QRb' or 'QRseq'.") 
    }
    if (is.list(QRobj) == FALSE){
        pandterm("QRobj is not a list") 
    }
    if (length(sapply(X=QRobj[[1]],FUN=length))==1){
        nqr <- 1 
    } else {
        nqr <- length(sapply(X=QRobj,FUN=length))
    }
    if (nqr>1){
        for (i in 1:nqr){
            if(is.null(QRobj[[i]]$method)){
                pandterm("The primary lists of QRobj does not contain 'method'") 
            }
            if(is.null(QRobj[[i]]$p)){
                pandterm("The primary lists of QRobj does not contain the quantile 'p'") 
            }
            if(is.null(QRobj[[i]]$betadraw)){
                pandterm("The primary lists of QRobj does not contain the beta estimates 'betadraw'") 
            }
        }
    } else {
        if(is.null(QRobj$method)){
            pandterm("QRobj does not contain 'method'") 
        }
        if(is.null(QRobj$p)){
            pandterm("QRobj does not contain the quantile 'p'") 
        }
        if(is.null(QRobj$betadraw)){
            pandterm("QRobj does not contain the beta estimates 'betadraw'") 
        }
    }


    # Calculate Bayes estimate and posterior credible interval if traceplot=FALSE 
    if (!traceplot){
	# In case of QRseq object do this
        if (nqr>1){
            k <- ncol(QRobj[[1]]$betadraw)
            ndraw <- nrow(QRobj[[1]]$betadraw)
            if(burnin>ndraw){
                pandterm("Burnin cannot be larger than number of mcmc draws R") 
            }
            res <- matrix(NA, nrow=nqr*k, ncol=3)
            rownames <- NULL
            colnames <- c("Bayes Estimate", as.character(credint)) 
            tmp <- rep(NA,3)
            iii <- 0
            for (i in 1:nqr){
                for (ii in 1:k){
                rownames <- c(rownames, paste("Quantile: ", QRobj[[i]]$p, "  Beta: ", ii, "   ", sep=""))
                tmp[1] <- mean(QRobj[[i]]$betadraw[burnin:ndraw,ii])
                tmp[2:3] <- quantile(QRobj[[i]]$betadraw[burnin:ndraw,ii], probs=credint)
                iii <- iii+1
                res[iii,] <- round(tmp, 3)
                }
            }
            colnames(res) <- colnames
            rownames(res) <- rownames
    
	# In case of QRb, QRc of QRc.AL object do this
        } else {
            k <- ncol(QRobj$betadraw)
            ndraw <- nrow(QRobj$betadraw)
            if(burnin>ndraw){
                pandterm("Burnin cannot be larger than number of mcmc draws R") 
            }
            res <- matrix(NA, nrow=k, ncol=3)
            rownames <- NULL
            colnames <- c("Bayes Estimate", as.character(credint)) 
            tmp <- rep(NA,3)
            for (i in 1:k){
                rownames <- c(rownames, paste("Quantile: ", QRobj$p, "  Beta: ", i, "   ", sep=""))
                tmp[1] <- mean(QRobj$betadraw[burnin:ndraw,i])
                tmp[2:3] <- quantile(QRobj$betadraw[burnin:ndraw,i], probs=credint)
                res[i,] <- round(tmp, 3)
            }
            colnames(res) <- colnames
            rownames(res) <- rownames
        }
    
        return(res)

    # Draw traceplots for all parameters 
    } else {
        # In case of QRseq object do this 
        if (nqr>1){
            k <- ncol(QRobj[[1]]$betadraw)
            ndraw <- nrow(QRobj[[1]]$betadraw)
            if(burnin>ndraw){
                pandterm("Burnin cannot be larger than number of mcmc draws R") 
            }
            for (i in 1:nqr){
                for (ii in 1:k){
                    plot(QRobj[[i]]$betadraw[burnin:ndraw,ii], typ="l", xlab="iteration", ylab="beta",
    		         main=paste("Quantile: ", QRobj[[i]]$p, " - Beta ", ii))
                    if ((i == nqr) & (ii == k)) break 
                    ans <- readline("Do you want to see the next traceplot (type 'yes' or 'no'):\n")
                    while ((ans != "yes") & (ans != "no")) ans <- readline("Incorrect input, type 'yes' or 'no':\n")
                    if (ans == "no") break 
                }
                if (ans == "no") break 
            }
	# In case of QRb, QRc of QRc.AL object do this
        } else {
            k <- ncol(QRobj$betadraw)
            ndraw <- nrow(QRobj$betadraw)
            if(burnin>ndraw){
                pandterm("Burnin cannot be larger than number of mcmc draws R") 
            }
            for (i in 1:k){
                plot(QRobj$betadraw[burnin:ndraw,i], typ="l", xlab="iteration", ylab="beta",
		     main=paste("Quantile: ", QRobj$p, " - Beta ", i))
                if (i == k) break 
                ans <- readline("Do you want to see the next traceplot (type 'yes' or 'no'):\n")
                while ((ans != "yes") & (ans != "no")) ans <- readline("Incorrect input, type 'yes' or 'no':\n")
                if (ans == "no") break 
            }
	}
    }
} 
