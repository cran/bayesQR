QRplot <- function(QRseq.obj, var=1, burnin=1, credint=c(0.025,0.975), 
                   main="", xlab="quantile",ylab="beta",xlim=c(0,1)){

  # Error handling
    pandterm = function(message) {
        stop(message, call. = FALSE)
    }
  	if (missing(QRseq.obj)){
	      pandterm("Requires QRseq.obj as provided by the QRseq function") 
	  }
  	if (is.list(QRseq.obj) == FALSE){
	      pandterm("QRseq.obj is not a list") 
	  }
		nqr <- length(sapply(X=QRseq.obj,FUN=length))
  	if (nqr < 2){
	      pandterm("QRseq.obj should contain the results of at least 2 estimated quantiles") 
	  }
		for (i in 1:nqr){
		    if(is.null(QRseq.obj[[i]]$method)){
	          pandterm("The primary lists of QRseq.obj do not contain 'method'") 
				}
		    if(is.null(QRseq.obj[[i]]$p)){
	          pandterm("The primary lists of QRseq.obj do not contain the quantile 'p'") 
				}
		    if(is.null(QRseq.obj[[i]]$betadraw)){
	          pandterm("The primary lists of QRseq.obj do not contain the beta estimates 'betadraw'") 
				}
		}


		# Create the quantile plot data
		plotdata <- matrix(NA,nrow=nqr,ncol=4)
		for (i in 1:nqr){
      ndraw <- nrow(QRseq.obj[[i]]$betadraw)
      if (burnin > ndraw) pandterm("Burnin cannot be larger than number of mcmc draws R") 
			# Estimated quantile (x-axis)
	    plotdata[i,1] <- QRseq.obj[[i]]$p  
			# Bayes estimate
	    plotdata[i,2] <- mean(QRseq.obj[[i]]$betadraw[burnin:ndraw,var])
			# Credible interval
	    plotdata[i,3:4] <- quantile(QRseq.obj[[i]]$betadraw[burnin:ndraw,var],probs=credint)
		}

    # plot axes/box in correct scale
		plot(x=NULL, y=NULL, xlim=xlim, ylim=c(min(plotdata[,2:4]),max(plotdata[,2:4])),
		     main=main, xlab=xlab, ylab=ylab)

    # calculate a small value (outside/below the box region)
    small <- min(plotdata[,2:4])-(max(plotdata[,2:4])-min(plotdata[,2:4]))

    # plot credible interval
		polygon(x=c(plotdata[1,1],plotdata[,1],plotdata[nqr,1]), 
		        y=c(small,plotdata[,4],small),col="grey",border=FALSE)

		polygon(x=c(plotdata[1,1],plotdata[,1],plotdata[nqr,1]), 
		        y=c(small,plotdata[,3],small),col="white",border=FALSE)

    # plot Bayes estimate
		points(x=plotdata[,1],y=plotdata[,2],typ="o",lty=2)

		# some aesthetic additions
		points(x=plotdata[,1],y=plotdata[,3],typ="l",col="darkgrey")
		points(x=plotdata[,1],y=plotdata[,4],typ="l",col="darkgrey")
		box(lwd=1.3,col="white")
		box(lwd=1.3,col="black")
}
