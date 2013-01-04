QRseq <- function(Method, Data, Prior, Mcmc){

  # Error handling
    pandterm = function(message) {
        stop(message, call. = FALSE)
    }
    if (missing(Method)) {
        pandterm("Requires Method argument -- should be 'QRc', 'QRc.AL', 'QRb' or 'QRb.AL'")
    }
    if ((Method!="QRb") & (Method!="QRc") & (Method!="QRb.AL") & (Method!="QRc.AL")) {
        pandterm("Method argument should be 'QRc', 'QRc.AL', 'QRb' or 'QRb.AL'")
    }
    if (missing(Data)) {
        pandterm("Requires Data argument -- list of y, X and p")
    }
    if (is.null(Data$p)) {
        pandterm("Requires Data element p")
    }
    pvec = Data$p
	nqr = length(pvec)
    if (nqr < 2) {
        pandterm("At least 2 quantiles should be estimated")
    }
    for (i in 1:nqr){
        if ((pvec[i]<=0)|(pvec[i]>=1)) {
            pandterm(paste("Unvalid quantile p[",i,"] =",pvec[i],sep=""))
        }
    }


  # define an empty object (list)
  out <- NULL 

  # estimate the models
  for (i in 1:nqr){
      # print information to console
          cat("************************************************","\n")
          cat("* Start estimating quantile ", i," of ", nqr, "in total *", "\n")
          cat("************************************************","\n")
      # set correct quantile and estimate model
	  Data$p <- pvec[i]
      if (Method=="QRc"){
          out[[i]] <- QRc(Data=Data,Prior=Prior,Mcmc=Mcmc)
      } else if (Method == "QRb"){
          out[[i]] <- QRb(Data=Data,Prior=Prior,Mcmc=Mcmc)
      } else if (Method == "QRc.AL"){
          out[[i]] <- QRc.AL(Data=Data,Prior=Prior,Mcmc=Mcmc)
      } else if (Method == "QRb.AL"){
          out[[i]] <- QRb.AL(Data=Data,Prior=Prior,Mcmc=Mcmc)
      }
  }
  return(out)
}
