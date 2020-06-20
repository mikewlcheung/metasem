calEffSizes <- function(model, n, Cov, Mean, lavaan.output=FALSE, ...) {

    if (missing (Mean)) {
        fit <- lavaan::sem(model, sample.cov=Cov, sample.nobs=n, sample.cov.rescale=FALSE, ...)
    } else {
        fit <- lavaan::sem(model, sample.cov=Cov, sample.mean=Mean, sample.nobs=n, sample.cov.rescale=FALSE, ...)
    }
    
    if (lavaan.output==FALSE) {
       
        ## Obtain the free parameters in the model
        x <- fit@Fit@x

        ## Obtain the sampling covariance matrix of the parameter estimates
        VCOV <- lavaan::vcov(fit)

        ## Get the effect sizes
        ES <- fit@Model@def.function(.x.=x)

        ## Compute the jacobian for 'defined parameters'
        JAC <- .lavJacobianD(func=fit@Model@def.function, x=x)

        ## Compute the sampling covariance matrix using delta method
        ES.VCOV <- JAC %*% VCOV %*% t(JAC)

        ## Add the variable names for ease of reference
        dimnames(ES.VCOV) <- list(names(ES), names(ES))

        fit <- list(ES=ES, VCOV=ES.VCOV)
    }
    fit
}


