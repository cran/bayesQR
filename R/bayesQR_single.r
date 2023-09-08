bayesQR_single <- function(formula=NULL, data=NULL, quantile=0.5, alasso=FALSE, normal.approx=NULL, ndraw=NULL, keep=1, prior=NULL){

	# Create function for error message
	pandterm <- function(message) {
	  stop(message, call. = FALSE)
	}

	# Extract input data
	if (is.null(formula)){
		pandterm("Formula has to be specified")
	}

	mf <- model.frame(formula=formula, data=data)
	X <- model.matrix(attr(mf,"terms"), data=mf)
	y <- model.response(mf)
	n <- length(y)
	nvar <- ncol(X)
	names <- colnames(X)

	# Check specified quantile
	if ((quantile<=0)|(quantile>=1)) {
		pandterm("Quantiles should be between zero and one")
	}

	# Check specified number of mcmc draws
	if (is.null(ndraw)) {
		pandterm("Number of mcmc draws (ndraw) has to be specified")
	}

	# Check specified number of mcmc draws
	if (keep>ndraw) {
		pandterm("'keep' cannot be larger than 'ndraw'")
	}

	# Determine method based on depvar
	dv <- unique(y)

	## Reset option flags
	QRc <- FALSE; QRc.AL <- FALSE; QRb <- FALSE; QRb.AL <- FALSE;

	## Determine method
	if (length(dv)>2){
		if (alasso) {
			QRc.AL <- TRUE
		} else {
			QRc <- TRUE
		}
	} else if (length(dv)==2){
		if (!identical(as.numeric(sort(dv)),c(0,1))){
			pandterm("The dependent variable should be coded as vector of zeros and ones")
		}
		if (alasso) {
			QRb.AL <- TRUE
		} else {
			QRb <- TRUE
		}
	} else {
		pandterm("Something is wrong with the dependent variable")
	}

	# QRc 
	#================================
	if (QRc) {
		# If normal approx was left unspecified, switch it on
		if (is.null(normal.approx)) normal.approx <- TRUE

		# If missing, set default prior 
		if (is.null(prior)) {
			prior <- prior(formula=formula, data=data, alasso=alasso)
		# Else check provided prior
		} else if (is(prior)!="bayesQR.prior"){
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
		} else if (prior$method!="QRc") {
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
		} else if (normal.approx & ((prior$shape0 != .01) | (prior$scale0 != .01))){
			warning("When normal.approx equals TRUE, then sigma is fixed and a prior for sigma is thus irrelevant")
		}

		# Assign correct variable types
		n <- as.integer(n)
		nvar <- as.integer(nvar)
		ndraw <- as.integer(ndraw)
		keep <- as.integer(keep)
		y <- as.double(y)
		quantile <- as.double(quantile)
		X <- as.double(X)
		beta0 <- as.double(prior$beta0)
		V0i <- as.double(chol2inv(chol(prior$V0)))
		shape0 <- as.double(prior$shape0)
		scale0 <- as.double(prior$scale0)
		normal <- as.logical(normal.approx)
		betadraw <- double(nvar*ndraw/keep)
		sigmadraw <- double(ndraw/keep)

		# Call Fortran routine
		fn_val <- .Fortran("QRc_mcmc", n, nvar, ndraw, keep, y, quantile, X, beta0, V0i, shape0, scale0, normal, betadraw, sigmadraw)

		# Normal approximation of posterior
		if (normal.approx){
			sigma.normal <- matrix(fn_val[[13]],nrow=ndraw/keep,ncol=nvar)
			sigma.normal <- sweep(sigma.normal,2,colMeans(sigma.normal))
			sigma.normal <- t(sigma.normal)%*%sigma.normal/(ndraw/keep)
			sigma.normal <- n/sqrt(n)*quantile*(1-quantile)*sigma.normal%*%(t(matrix(X,nrow=n))%*%matrix(X,nrow=n))%*%sigma.normal
			sigma.normal <- as.vector(sigma.normal)
		}

		# Return bayesQR object
		out <- list(method="QRc",
			    normal.approx=normal.approx,
			    quantile=quantile,
			    names=names,
			    betadraw=matrix(fn_val[[13]],nrow=ndraw/keep,ncol=nvar),
			    sigmadraw=fn_val[[14]],
			    sigma.normal=ifelse(rep(normal.approx,nvar*nvar),sigma.normal,NA))

	# QRc.AL
	#================================
	} else if (QRc.AL) {
		# If normal approx was left unspecified, switch it on
		if (is.null(normal.approx)) normal.approx <- TRUE

		# If missing, set default prior 
		if (is.null(prior)) {
			prior <- prior(formula=formula, data=data, alasso=alasso)
		# Else check provided prior
		} else if (is(prior)!="bayesQR.prior"){
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
		} else if (prior$method!="QRc.AL") {
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
		} else if (normal.approx & ((prior$shape0 != .01) | (prior$scale0 != .01))){
			warning("when normal.approx equals TRUE, then sigma is fixed and a prior for sigma is thus irrelevant")
		}

		# Assign correct variable types
		n <- as.integer(n)
		nvar <- as.integer(nvar)
		ndraw <- as.integer(ndraw)
		keep <- as.integer(keep)
		y <- as.double(y)
		quantile <- as.double(quantile)
		x <- as.double(X)
		a <- as.double(prior$a)
		b <- as.double(prior$b)
		c <- as.double(prior$c)
		d <- as.double(prior$d)
		normal <- as.logical(normal.approx)
		betadraw <- double(nvar*ndraw/keep)
		sigmadraw <- double(ndraw/keep)

		# Call Fortran routine
		fn_val <- .Fortran("QRc_AL_mcmc",n, nvar, ndraw, keep, y, quantile, x, a, b, c, d, normal, betadraw, sigmadraw) 

		# Normal approximation of posterior
		if (normal.approx){
			sigma.normal <- matrix(fn_val[[13]],nrow=ndraw/keep,ncol=nvar)
			sigma.normal <- sweep(sigma.normal,2,colMeans(sigma.normal))
			sigma.normal <- t(sigma.normal)%*%sigma.normal/(ndraw/keep)
			sigma.normal <- n/sqrt(n)*quantile*(1-quantile)*sigma.normal%*%(t(matrix(X,nrow=n))%*%matrix(X,nrow=n))%*%sigma.normal
			sigma.normal <- as.vector(sigma.normal)
		}

		# Return bayesQR object
		out <- list(method="QRc.AL",
			    normal.approx=normal.approx,
			    quantile=quantile,
			    names=names,
			    betadraw=matrix(fn_val[[13]],nrow=ndraw/keep, ncol=nvar),
			    sigmadraw=fn_val[[14]],
			    sigma.normal=ifelse(rep(normal.approx,nvar*nvar),sigma.normal,NA))

	# QRb
	#================================
	} else if (QRb) {
		# If normal approx was left unspecified, switch it off
		if (is.null(normal.approx)) normal.approx <- FALSE 

		# If missing, set default prior 
		if (is.null(prior)) {
			prior <- prior(formula=formula, data=data, alasso=alasso)
		# Else check provided prior
		} else if (is(prior)!="bayesQR.prior"){
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
		} else if (prior$method!="QRb") {
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
		}
		
		# If user requested the normal approximation, give warning
		if (normal.approx){
			warning("Normal approximation of the posterior is not supported for binary dependent variables")
		}

		# Assign correct variable types
		n <- as.integer(n)
		nvar <- as.integer(nvar)
		ndraw <- as.integer(ndraw)
		keep <- as.integer(keep)
		y <- as.integer(y)
		quantile <- as.double(quantile)
		X <- as.double(X)
		beta0 <- as.double(prior$beta0)
		V0i <- as.double(chol2inv(chol(prior$V0)))
		betadraw <- double(nvar*ndraw/keep)
		sigmadraw <- double(ndraw/keep)
		
		# Call Fortran routine
		fn_val <- .Fortran("QRb_mcmc", n, nvar, ndraw, keep, y, quantile, X, beta0, V0i, betadraw)

		# Return bayesQR object
		out <- list(method="QRb",
			    normal.approx=FALSE,
			    quantile=quantile, 
			    names=names,
			    betadraw=matrix(fn_val[[10]], nrow=ndraw/keep,ncol=nvar),
			    sigma.normal=rep(NA,nvar*nvar)
			    )
	
	# QRb.AL
	#================================
	} else if (QRb.AL) {

		# If normal approx was left unspecified, switch it off
		if (is.null(normal.approx)) normal.approx <- FALSE 

		# If missing, set default prior 
		if (is.null(prior)) {
			prior <- prior(formula=formula, data=data, alasso=alasso)
		# Else check provided prior
		} else if (is(prior)!="bayesQR.prior"){
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
		} else if (prior$method!="QRb.AL") {
			pandterm("Incorrect prior: type '?prior' for more information about how to define priors")
		} 

		# If user requested the normal approximation, give warning
		if (normal.approx){
			warning("Normal approximation of the posterior is not supported for binary dependent variables")
		}
		
		# Assign correct variable types
		n <- as.integer(n)
		nvar <- as.integer(nvar)
		ndraw <- as.integer(ndraw)
		keep <- as.integer(keep)
		y <- as.integer(y)
		quantile <- as.double(quantile)
		x <- as.double(X)
		c <- as.double(prior$c)
		d <- as.double(prior$d)
		betadraw <- double(nvar*ndraw/keep)
		
		# Call Fortran routine
		fn_val <- .Fortran("QRb_AL_mcmc",n, nvar, ndraw, keep, y, quantile, x, c, d, betadraw) 

		# Return bayesQR object
		out <- list(method="QRb.AL",
			    normal.approx=FALSE,
			    quantile=quantile,
			    names=names,
			    betadraw=matrix(fn_val[[10]], nrow=ndraw/keep, ncol=nvar),
			    sigma.normal=rep(NA,nvar*nvar)
			    )
		
	} else {
		pandterm("Something is wrong with the dependent variable")
	}
return(out)
} 
