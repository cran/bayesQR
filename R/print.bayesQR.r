print.bayesQR <- function(x, digits=3, ...){

	# Create summary object of x
	x <- summary(x)

	# Number of quantile regression summaries to print
	nqr <- length(x)

	# Loop trough every quantile regression 
	for (i in 1:nqr){
		QRsub <- x[[i]]
		cat("\n")
		if (QRsub$method=="QRc"){
			cat("Type of dependent variable: continuous\n")
			cat("Lasso variable selection: no\n")
		} else if (QRsub$method=="QRc.AL"){
			cat("Type of dependent variable: continuous\n")
			cat("Lasso variable selection: yes\n")
		} else if (QRsub$method=="QRb"){
			cat("Type of dependent variable: binary\n")
			cat("Lasso variable selection: no\n")
		} else if (QRsub$method=="QRb.AL"){
			cat("Type of dependent variable: binary\n")
			cat("Lasso variable selection: yes\n")
		}
		if (QRsub$normal.approx){
			cat("Normal approximation of posterior: yes\n")
		} else {
			cat("Normal approximation of posterior: no\n")
		}
		cat(paste("Estimated quantile: ",QRsub$quantile),"\n")
		cat("\n")
		cat("\n")
		cat("Coefficients (Bayes estimates):\n")
		cat("\n")
		print(QRsub$betadraw[,1],digits=digits)
		if ((nqr>1)&(i<nqr)) cat("*****************************************\n")
	}
}
