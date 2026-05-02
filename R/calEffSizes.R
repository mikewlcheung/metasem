#' Calculate Effect Sizes using lavaan Models
#'
#' It calculates effect sizes with Delta Method by formulating the effect sizes
#' as functions of SEM in lavaan.
#'
#'
#' @param model A lavaan model. Effect sizes are defined as functions of SEM
#' parameters with \code{:=}.
#' @param data A data frame of the observed variables. If it is \code{NULL},
#' summary statistics are required.
#' @param n Sample sizes
#' @param Cov A covariance matrix or a list of covariance matrices.
#' @param Mean Optional sample means.
#' @param group A character of the variable name in the data frame defining the
#' groups in a multiple group analysis.
#' @param lavaan.output If \code{TRUE}, it returns the fitted object instead of
#' the effect sizes and their sampling covariance matrix.
#' @param warn If \code{FALSE}, it suppresses lavaan related warnings.
#' @param \dots Further arguments passed to \code{\link[lavaan]{sem}}.
#' @return Effect sizes and their sampling covariance matrix or a lavaan fitted
#' object.
#' @note The input matrices are treated as covariance matrices unless there are
#' explicit constraints in the model.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{smdMES}}, \code{\link[metaSEM]{smdMTS}}
#' @references Cheung, M. W.-L. (2015). \emph{Meta-analysis: A structural
#' equation modeling approach}. Chichester, West Sussex: John Wiley & Sons,
#' Inc.
#'
#' Cheung, M. W.-L. (2018). Computing multivariate effect sizes and their
#' sampling covariance matrices with structural equation modeling: Theory,
#' examples, and computer simulations. \emph{Frontiers in Psychology},
#' \bold{9}(1387). https://doi.org/10.3389/fpsyg.2018.01387
#' @keywords meta-analysis
#' @examples
#'
#' \donttest{ 
#' ## Select ATT, Bi, and BEH 
#' obs.vars <- c("BEH", "BI", "ATT")
#'
#' ## Select one study from Cooke16 for illustration
#' my.cor <- Cooke16$data[[4]][obs.vars, obs.vars]
#' my.n  <- Cooke16$n[4]
#'
#' ## Effect sizes: indirect effect and direct effect
#' model <- "BEH ~ c*ATT + b*BI
#'           BI ~ a*ATT
#'           ## Indirect effect
#'           Ind := a*b
#'           Dir := c"
#'
#' calEffSizes(model=model, n=my.n, Cov=my.cor, lavaan.output=FALSE)
#'
#' ## Return the lavaan fitted model
#' fit <- calEffSizes(model=model, n=my.n, Cov=my.cor, lavaan.output=TRUE)
#' lavaan::summary(fit)
#'
#' lavaan::parameterestimates(fit)  
#' }
#'
calEffSizes <- function(model, data=NULL, n, Cov, Mean=NULL, group=NULL,
                        lavaan.output=FALSE, warn=FALSE, ...) {

    ## When raw data are present
    if (!is.null(data)) {
        fit <- lavaan::sem(model, data=data, group=group, warn=warn, ...)
    } else {
        ## Summary statistics as inputs     
        if (is.null (Mean)) {
          fit <- lavaan::sem(model, sample.cov=Cov, sample.nobs=n, group=group,
                             sample.cov.rescale=FALSE, warn=warn, ...)
        } else {
          fit <- lavaan::sem(model, sample.cov=Cov, sample.mean=Mean,
                             sample.nobs=n, group=group,
                             sample.cov.rescale=FALSE, warn=warn, ...)
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


