calEffSizes <- function(model, data=NULL, n, Cov, Mean=NULL, group=NULL, lavaan.output=FALSE, ...) {

    ## When raw data are present
    if (!is.null(data)) {
        fit <- lavaan::sem(model, data=data, group=group, ...)
    } else {
        ## Summary statistics as inputs     
        if (is.null (Mean)) {
            fit <- lavaan::sem(model, sample.cov=Cov, sample.nobs=n, group=group, sample.cov.rescale=FALSE, ...)
        } else {
            fit <- lavaan::sem(model, sample.cov=Cov, sample.mean=Mean, sample.nobs=n, group=group, sample.cov.rescale=FALSE, ...)
        }
    }
    
    if (lavaan.output==FALSE) {
       
        ## Get the free parameters in the model
        x <- fit@Fit@x

        ## Get the sampling covariance matrix of the parameter estimates
        VCOV <- lavaan::vcov(fit)

        ## Compute the effect sizes
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


