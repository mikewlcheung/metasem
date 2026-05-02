## This function is directly copied from lavaan to avoid a warning from calling lavaan:::lavJacobianD
.lavJacobianD <- function (func, x, h = sqrt(.Machine$double.eps), ...) 
{
    f0 <- func(x, ...)
    nres <- length(f0)
    nvar <- length(x)
    h <- pmax(h, abs(h * x))
    tmp <- x + h
    h <- (tmp - x)
    dx <- matrix(as.numeric(NA), nres, nvar)
    for (p in seq_len(nvar)) {
        dx[, p] <- (func(x + h * (seq.int(nvar) == p), ...) - 
            func(x, ...))/h[p]
    }
    dx
}



#' Compute Effect Sizes for Multiple Treatment Studies
#' 
#' It computes the standardized mean differences and their asymptotic sampling
#' covariance matrix for \emph{k} multiple treatment studies. The first group
#' is assumed as the control group.
#' 
#' Gleser and Olkin (2009) introduce formulas to calculate the standardized
#' mean differences and their sampling covariance matrix for multiple treatment
#' studies under the assumption of homogeneity of the covariance matrix. This
#' function uses a structural equation modeling (SEM) approach introduced in
#' Chapter 3 of Cheung (2015) to calculate the same estimates. The SEM approach
#' is more flexible in three ways: (1) it allows homogeneity of variances or
#' not; (2) it allows users to test the assumption of homogeneity of variances
#' by checking the fitted \code{\link[lavaan]{lavaan-class}} object; and (3) it
#' may calculate all pairwise comparisons.
#' 
#' @param m A vector of \emph{k} sample means.
#' @param v A vector of \emph{k} sample variances.
#' @param n A vector of \emph{k} sample sizes.
#' @param homogeneity If it is \code{variance} (the default), homogeneity of
#' variances is assumed. The common standard deviation is used as the
#' standardizer in calculating the effect sizes. If it is \code{none},
#' homogeneity of variances is not assumed. The standard deviation of the first
#' group is used as the standardizer in calculating the effect sizes.
#' @param bias.adjust If it is \code{TRUE} (the default), the effect sizes are
#' adjusted for small bias by multiplying \eqn{1-3/(4*(n1+n2)-9)}.
#' @param all.comparisons If it is \code{FALSE} (the default), all groups
#' (except the first group) are compared against the first group. If it is
#' \code{TRUE}, all pairwise comparisons are calculated. This may be useful in
#' network meta-analysis.
#' @param list.output If it is \code{TRUE} (the default), the effect sizes and
#' their sampling covariance matrix are outputed as a list. If it is
#' \code{FALSE}, they will be stacked into a vector.
#' @param lavaan.output If it is \code{FALSE} (the default), the effect sizes
#' and its sampling covariance matrix are reported. If it is \code{TRUE}, it
#' outputs the fitted \code{\link[lavaan]{lavaan-class}} object.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{Gleser94}}, \code{\link[metaSEM]{smdMES}},
#' \code{\link[metaSEM]{calEffSizes}}
#' @references Cheung, M. W.-L. (2015). \emph{Meta-analysis: A structural
#' equation modeling approach}. Chichester, West Sussex: John Wiley & Sons,
#' Inc.
#' 
#' Cheung, M. W.-L. (2018). Computing multivariate effect sizes and their
#' sampling covariance matrices with structural equation modeling: Theory,
#' examples, and computer simulations. \emph{Frontiers in Psychology},
#' \bold{9}(1387). https://doi.org/10.3389/fpsyg.2018.01387
#' 
#' Gleser, L. J., & Olkin, I. (2009). Stochastically dependent effect sizes. In
#' H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of
#' research synthesis and meta-analysis}. (2nd ed., pp. 357-376). New York:
#' Russell Sage Foundation.
#' @keywords meta-analysis
#' @examples
#' 
#' \donttest{  
#' ## Sample means for groups 1 to 3
#' m <- c(5,7,9)
#' 
#' ## Sample variances
#' v <- c(10,11,12)
#' 
#' ## Sample sizes
#' n <- c(50,52,53)
#' 
#' ## Assuming homogeneity of variances
#' smdMTS(m, v, n, homogeneity = "var", bias.adjust=TRUE, all.comparisons=FALSE,
#'        lavaan.output=FALSE)
#' 
#' ## Not assuming homogeneity of variances and comparing all pairwise groups
#' ## Please note that the SD of the first group is used as the standardizer    
#' smdMTS(m, v, n, homogeneity = "none", bias.adjust=TRUE, all.comparisons=TRUE,
#'        lavaan.output=FALSE)
#' 
#' ## Output the fitted lavaan model
#' ## It provides a likelihood ratio test to test the null hypothesis of
#' ## homogeneity of variances.    
#' fit <- smdMTS(m, v, n, homogeneity = "var", bias.adjust=FALSE, all.comparisons=FALSE,
#'               lavaan.output=TRUE)
#' 
#' lavaan::summary(fit)
#'     
#' lavaan::parameterestimates(fit)
#' }
#' 
smdMTS <- function(m, v, n, homogeneity=c("variance", "none"), bias.adjust=TRUE, 
                   all.comparisons=FALSE, list.output=TRUE, lavaan.output=FALSE) {

    ## if (!requireNamespace("lavaan", quietly=TRUE))
    ##   stop("\"lavaan\" package is required for this function.")

    if ( var(c(length(m), length(v), length(n))) !=0 ) {
        stop("The lengths of m, v, and n are not the same.\n")
    }

    ## ind: Complete data without NA
    ind <- m&n&v
    ## Only TRUE or FALSE, no NA
    ind[is.na(ind)] <- FALSE
    k_tot <- length(ind)

    if (all.comparisons) {
        ind.y <- names.y <- rep(NA, k_tot*(k_tot-1)/2)
        p <- 1
        for (i in 1:(k_tot-1))
            for (j in (i+1):k_tot) {
                ind.y[p] <- ind[i]&ind[j]
                names.y[p] <- paste0("y",j,"_",i)
                p <- p+1
        }
    } else {
        ind.y <- names.y <- rep(NA, k_tot-1)
        for (i in seq_len(k_tot-1)) {
            ind.y[i] <- ind[1]&ind[i+1]
            names.y[i] <- paste0("y",i+1,"_1")
        }
    }

    ## Remove NA
    m <- m[ind]
    v <- v[ind]
    n <- n[ind]

    ## No. of groups with complete data
    k <- length(m)

    homogeneity <- match.arg(homogeneity)
    if (k<2) stop("Number of groups must be at least 2.\n")
    if (k>2 & homogeneity=="none" & all.comparisons==TRUE) {
    ## The first group is the control group.
        warning("The standard deviation of the first group is used to calculate the effect sizes for all comparisons.\n")
    }

    ## A list of variance matrices
    Var <- lapply(v, function(x) matrix(x, dimnames = list("x", "x")))
    ## A list of mean vectors
    Mean <- lapply(m, function(x) {names(x) <- "x"; x})

    ## Assuming homogeneity of variances by using the same label "s1" for standard deviation
    if (homogeneity=="variance") {
        model1 <- paste0("lat =~ c(", paste0(rep("s1", k), collapse=","), 
                         ")*x+start(",paste0(rep(sqrt(mean(v)), length(v)), collapse=","),")*x\n")
    } else {
        model1 <- paste0("lat =~ c(", paste0("s", seq_len(k), collapse=","),
                         ")*x+start(",paste0(sqrt(v),collapse=","),")*x\n")
    }

    ## Means: m1, m2, m3
    model2 <- paste0("x ~ c(", paste0("m", seq_len(k), collapse=","), 
                     ")*1 + start(", paste0(m,collapse = ","), ")*1\n")

    ## Bias adjustment factor
    cm <- function(n1, n2) 1-3/(4*(n1+n2)-9)

    ## Functions of parameters
    model3 <- list()
    index <- 1
    if (all.comparisons) {
        for (i in 1:(k-1))
            for (j in (i+1):k) {
                if (bias.adjust==TRUE) {
                    model3[index] <- paste0("y",j,"_",i, " := ", cm(n[i],n[j]),"*(m",j,"-m",i,")/s1\n")
                } else {
                    model3[index] <- paste0("y",j,"_",i, " := (m",j,"-m",i,")/s1\n")
                }
                index <- index+1
            }
    } else {
        for (i in 1:(k-1))
            if (bias.adjust==TRUE) {
                model3[i] <- paste0("y",i+1,"_1 := ", cm(n[1],n[i+1]), "*(m", i+1, "-m1)/s1\n")
            } else {
                model3[i] <- paste0("y",i+1,"_1 := (m", i+1, "-m1)/s1\n")
            }
    }

    model <- paste0(model1, model2, do.call(paste0, model3))

    fit <- lavaan::sem(model, sample.cov=Var, sample.mean=Mean, std.lv=TRUE,
                       sample.nobs=n, sample.cov.rescale=FALSE)

    if (lavaan.output) {
        out <- fit
    } else {
        ## Obtain the free parameters in the model
        x <- fit@Fit@x
        ## Compute the multiple effect sizes
        y.complete <- fit@Model@def.function(.x.=x)
        ## Compute the jacobian for the 'defined parameters'
        JAC <- .lavJacobianD(func=fit@Model@def.function, x=x)
        ## Compute the sampling covariance matrix using delta method
        V.complete <- JAC %*% lavaan::inspect(fit, what="vcov") %*% t(JAC)

        ## Add the NA
        y <- rep(NA, length(ind.y))
        V <- matrix(NA, ncol=length(ind.y), nrow=length(ind.y))
        y[ind.y] <- y.complete
        V[ind.y, ind.y] <- V.complete

        ## Add the variable names for ease of reference
        names(y) <- names.y
        dimnames(V) <- list(names.y, names.y)
        
        if (list.output) {
            out <- list(y=y, V=V)
        } else {
            v <- vech(V)
            names.v <- paste("C(", outer(names.y, names.y, paste, sep = " "), ")", sep="")
            names(v) <- vech(matrix(names.v, ncol=length(names.y), nrow=length(names.y)))
            out <- c(y,v)
        }        
    }
  
    out
}




#' Compute Effect Sizes for Multiple End-point Studies
#' 
#' It computes the standardized mean differences and their asymptotic sampling
#' covariance matrix for two multiple end-point studies with \emph{p} effect
#' sizes.
#' 
#' Gleser and Olkin (2009) introduce formulas to calculate the standardized
#' mean differences and their sampling covariance matrix for multiple end-point
#' studies under the assumption of homogeneity of the covariance matrix. This
#' function uses a structural equation modeling (SEM) approach introduced in
#' Chapter 3 of Cheung (2015) to calculate the same estimates. The SEM approach
#' is more flexible in two ways: (1) it allows homogeneity of covariance or
#' correlation matrices or not; and (2) it allows users to test this assumption
#' by checking the fitted \code{\link[lavaan]{lavaan-class}} object.
#' 
#' @param m1 A vector of \emph{p} sample means of the first group.
#' @param m2 A vector of \emph{p} sample means of the second group.
#' @param V1 A \emph{p} by \emph{p} sample covariance matrix of the first
#' group.
#' @param V2 A \emph{p} by \emph{p} sample covariance matrix of the second
#' group.
#' @param n1 The sample size of the first group.
#' @param n2 The sample size of the second group.
#' @param homogeneity If it is \code{covariance} (the default), homogeneity of
#' covariance matrices is assumed. The common standard deviations are used as
#' the standardizers in calculating the effect sizes. If it is
#' \code{correlation}, homogeneity of correlation is not assumed. The standard
#' deviations of the first group are used as the standardizer in calculating
#' the effect sizes. If it is \code{none}, no homogeneity assumption is made.
#' The standard deviations of the first group are used as the standardizer in
#' calculating the effect sizes.
#' @param bias.adjust If it is \code{TRUE} (the default), the effect sizes are
#' adjusted for small bias by multiplying \eqn{1-3/(4*(n1+n2)-9)}.
#' @param list.output If it is \code{TRUE} (the default), the effect sizes and
#' their sampling covariance matrix are outputed as a list. If it is
#' \code{FALSE}, they will be stacked into a vector.
#' @param lavaan.output If it is \code{FALSE} (the default), the effect sizes
#' and its sampling covariance matrix are reported. If it is \code{TRUE}, it
#' outputs the fitted \code{\link[lavaan]{lavaan-class}} object.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{Gleser94}}, \code{\link[metaSEM]{smdMTS}},
#' \code{\link[metaSEM]{calEffSizes}}
#' @references Cheung, M. W.-L. (2015). \emph{Meta-analysis: A structural
#' equation modeling approach}. Chichester, West Sussex: John Wiley & Sons,
#' Inc.
#' 
#' Cheung, M. W.-L. (2018). Computing multivariate effect sizes and their
#' sampling covariance matrices with structural equation modeling: Theory,
#' examples, and computer simulations. \emph{Frontiers in Psychology},
#' \bold{9}(1387). https://doi.org/10.3389/fpsyg.2018.01387
#' 
#' Gleser, L. J., & Olkin, I. (2009). Stochastically dependent effect sizes. In
#' H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of
#' research synthesis and meta-analysis}. (2nd ed., pp. 357-376). New York:
#' Russell Sage Foundation.
#' @keywords meta-analysis
#' @examples
#' 
#' \donttest{    
#' ## Sample means for the two constructs in Group 1  
#' m1 <- c(2.5, 4.5)
#' 
#' ## Sample means for the two constructs in Group 2     
#' m2 <- c(3, 5)
#' 
#' ## Sample covariance matrix in Group 1    
#' V1 <- matrix(c(3,2,2,3), ncol=2)
#' 
#' ## Sample covariance matrix in Group 2
#' V2 <- matrix(c(3.5,2.1,2.1,3.5), ncol=2)
#' 
#' ## Sample size in Group 1
#' n1 <- 20
#' 
#' ## Sample size in Group 2    
#' n2 <- 25
#' 
#' ## SMD with the assumption of homogeneity of covariance matrix    
#' smdMES(m1, m2, V1, V2, n1, n2, homogeneity="cov", bias.adjust=TRUE,
#'        lavaan.output=FALSE)
#' 
#' ## SMD with the assumption of homogeneity of correlation matrix 
#' smdMES(m1, m2, V1, V2, n1, n2, homogeneity="cor", bias.adjust=TRUE,
#'        lavaan.output=FALSE)
#' 
#' ## SMD without any assumption of homogeneity
#' smdMES(m1, m2, V1, V2, n1, n2, homogeneity="none", bias.adjust=TRUE,
#'        lavaan.output=FALSE)
#' 
#' ## Output the fitted lavaan model
#' ## It provides a likelihood ratio test to test the null hypothesis of
#' ## homogeneity of variances.     
#' fit <- smdMES(m1, m2, V1, V2, n1, n2, homogeneity="cov", bias.adjust=TRUE,
#'               lavaan.output=TRUE)
#' 
#' lavaan::summary(fit)
#' 
#' lavaan::parameterestimates(fit)
#' }
#' 
smdMES <- function(m1, m2, V1, V2, n1, n2, homogeneity=c("covariance", "correlation", "none"), 
                   bias.adjust=TRUE, list.output=TRUE, lavaan.output=FALSE) {

    ## if (!requireNamespace("lavaan", quietly=TRUE))
    ##   stop("\"lavaan\" package is required for this function.")

    if ( var(c(length(m1), length(m2), nrow(V1), ncol(V1), nrow(V2), ncol(V2))) !=0 ) {
        stop("Dimensions of the inputs are not the same.\n")
    }

    ## ind: Complete data without NA (only check m1 and m2)
    ind <- m1&m2
    ## Only TRUE or FALSE, no NA
    ind[is.na(ind)] <- FALSE
    names.y <- paste0("y", seq_along(ind))

    ## Remove NA
    m1 <- m1[ind]
    m2 <- m2[ind]
    V1 <- V1[ind, ind]
    V2 <- V2[ind, ind]

    ## No. of variables
    p <- length(m1)

    var.names <- paste0("x", seq_len(p))
    Var <-list(V1, V2)
    Var <- lapply(Var, function(x) {dimnames(x) <- list(var.names, var.names); x})
    Mean <- list(m1, m2)
    Mean <- lapply(Mean, function(x) {names(x) <- var.names; x})

    homogeneity <- match.arg(homogeneity)

    ## Bias adjustment factor
    cm <- function(n1, n2) 1-3/(4*(n1+n2)-9)

    model4 <- model3 <- model2 <- model1 <- list()

    for (i in seq_len(p)) {
        if ( homogeneity == "covariance" ) {
            ## Assuming homogeneity of variances by using the same label "s1"
            ## starting value for the average standard deviations
            ## equality on covariance is added in next block
            V <- (V1+V2)/2
            model1[i] <- paste0("lat",i," =~ c(s",i,"_1, s",i,"_1)*x",i,
                                "+start(",sqrt(V[i,i]), ",",sqrt(V[i,i]), ")*x",i,"\n")
        } else {
            model1[i] <- paste0("lat",i," =~ c(s",i,"_1, s",i,"_2)*x",i,
                                "+start(",sqrt(V1[i,i]),",",sqrt(V2[i,i]),")*x",i,"\n")
        }

        ## Means
        model2[i] <- paste0("x",i," ~ c(m",i,"_1, m",i,"_2)*1+start(", m1[i],",",m2[i],")*1\n")
        ## Error variances fixed at 0
        model3[i] <- paste0("x",i," ~~ 0*x",i,"\n")

        ## Effect sizes
        if (bias.adjust==TRUE) {
            model4[i] <- paste0("y",i," := ", cm(n1,n2), "*(m",i,"_2-m",i,"_1)/s",i,"_1\n")
        } else {
            model4[i] <- paste0("y",i," := (m",i,"_2-m",i,"_1)/s",i,"_1\n")
        }
    }

    model <- paste0(do.call(paste0, list(model1, model2, model3, model4)), collapse = "")

    ## Homogeneity of correlation or covariance matrices
    if ( homogeneity != "none" ) {
        fit <- lavaan::sem(model, sample.cov=Var, sample.mean=Mean, group.equal="lv.covariances",
                           std.lv=TRUE, sample.nobs=c(n1,n2), sample.cov.rescale=FALSE)
    } else {
        fit <- lavaan::sem(model, sample.cov=Var, sample.mean=Mean, 
                   std.lv=TRUE, sample.nobs=c(n1,n2), sample.cov.rescale=FALSE)
    }

    if (lavaan.output) {
        out <- fit
    } else {

        ## Obtain the free parameters in the model
        x <- fit@Fit@x
        ## Compute the multiple effect sizes
        y.complete <- fit@Model@def.function(.x.=x)
        ## Compute the jacobian for the 'defined parameters'
        JAC <- .lavJacobianD(func=fit@Model@def.function, x=x)
        ## Compute the sampling covariance matrix using delta method
        V.complete <- JAC %*% lavaan::inspect(fit, what="vcov") %*% t(JAC)

        ## Add the NA
        y <- rep(NA, length(ind))
        V <- matrix(NA, ncol=length(ind), nrow=length(ind))
        y[ind] <- y.complete
        V[ind, ind] <- V.complete

        ## Add the variable names for ease of reference
        names(y) <- names.y
        dimnames(V) <- list(names.y, names.y)

        if (list.output) {
            out <- list(y=y, V=V)
        } else {
            v <- vech(V)
            names.v <- paste("C(", outer(names.y, names.y, paste, sep = " "), ")", sep="")
            names(v) <- vech(matrix(names.v, ncol=length(names.y), nrow=length(names.y)))
            out <- c(y,v)
        }
    }
    
  out
}
