QRb.pred <- function(QRseq.obj,X,burnin=1){

  # Error handling
    pandterm = function(message) {
      stop(message, call. = FALSE)
    }
    if (missing(QRseq.obj)){
      pandterm("Requires QRseq.obj as provided by the QRseq function") 
    }
    if (is.list(QRseq.obj) == FALSE){
      pandterm("Requires QRseq.obj as provided by the QRseq function") 
    }
    nqr <- length(sapply(X=QRseq.obj,FUN=length))
    if (nqr < 9){
      pandterm("QRseq.obj should contain the results of at least 9 estimated quantiles") 
    }
    tmp <- 0
    for (i in 1:nqr){
      if(sapply(QRseq.obj,length)[i] != 3){
        pandterm("QRseq.obj has invalid dimensions") 
      }
      if( (QRseq.obj[[i]]$method != "QRb") & (QRseq.obj[[i]]$method != "QRb.AL") ){
        pandterm("QRseq.obj is not based on the 'QRb' or 'QRb.AL' method") 
      }
      if(ncol(QRseq.obj[[i]]$betadraw) != ncol(X)){
        pandterm("QRseq.obj has different number of estimated parameters than X has predictors") 
      }
      if(QRseq.obj[[i]]$p < tmp){
        pandterm("Lists in QRseq.obj must be sorted non-decreasingly p") 
      }
      tmp <- QRseq.obj[[i]]$p
    }

  # Make prediction of latent utility based on the Bayes estimate 
  n <- nrow(X)
  preds <- matrix(NA,nrow=n,ncol=nqr)
  for (i in 1:nqr){
    ndraw <- nrow(QRseq.obj[[i]]$betadraw)
    if (burnin > ndraw) pandterm("Burnin cannot be larger than number of mcmc draws R") 
    preds[,i] <- X%*%colMeans(QRseq.obj[[i]]$betadraw[burnin:ndraw,])
  }

  # Find interval that contains zero
  preds <- cbind(0,preds)
  preds <- t(apply(preds,FUN=sort,MARGIN=1))
  preds <- t(apply(preds,FUN="==",MARGIN=1,0))
  preds <- apply(preds,FUN=which,MARGIN=1)
  
  # Link utilities to probabilities
  pvec <- rep(NA,nqr)
  for (i in 1:nqr){
    pvec[i] <- QRseq.obj[[i]]$p
  }
  pvec <- (c(0,pvec)+c(pvec,1))/2

  # Pr(y=1) = 1 - p
  return(1-pvec[preds])
}
