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

smdMTS <- function(m, v, n, homogeneity=c("variance", "none"), bias.adjust=TRUE, 
                   all.comparisons=FALSE, lavaan.output=FALSE) {

  if (!requireNamespace("lavaan", quietly=TRUE))
    stop("\"lavaan\" package is required for this function.")
  
  if ( var(c(length(m), length(v), length(n))) !=0 ) {
    stop("The lengths of m, v, and n are not the same.\n")
  }
  
  ## No. of groups
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
  y <- fit@Model@def.function(.x.=x)
  ## Compute the jacobian for the 'defined parameters'
  JAC <- .lavJacobianD(func=fit@Model@def.function, x=x)
  ## Compute the sampling covariance matrix using delta method
  V <- JAC %*% lavaan::inspect(fit, what="vcov") %*% t(JAC)
  ## Add the variable names for ease of reference
  dimnames(V) <- list(names(y), names(y))
  out <- list(y=y, V=V)
  }
  
  out
}


smdMES <- function(m1, m2, V1, V2, n1, n2, homogeneity=c("covariance", "correlation", "none"), 
                   bias.adjust=TRUE, lavaan.output=FALSE) {

  if (!requireNamespace("lavaan", quietly=TRUE))
    stop("\"lavaan\" package is required for this function.")
  
  if ( var(c(length(m1), length(m2), nrow(V1), ncol(V1), nrow(V2), ncol(V2))) !=0 ) {
    stop("Dimensions of the inputs are not the same.\n")
  }
  
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
    
  ## Assuming homogeneity of variances by using the same label "s1"
  model4 <- model3 <- model2 <- model1 <- list()
    
  for (i in seq_len(p)) {        
    if ( homogeneity == "covariance" ) {
      ## starting value for the average standard deviations
      ## equality on covariance is added in next block
      V <- (V1+V2)/2
      model1[i] <- paste0("lat",i," =~ c(s",i,"_1, s",i,"_1)*x",i,"+start(",sqrt(V[i,i]), ",",sqrt(V[i,i]), ")*x",i,"\n")
    } else {
      model1[i] <- paste0("lat",i," =~ c(s",i,"_1, s",i,"_2)*x",i,"+start(",sqrt(V1[i,i]),",",sqrt(V2[i,i]),")*x",i,"\n")
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
  y <- fit@Model@def.function(.x.=x)
  ## Compute the jacobian for the 'defined parameters'
  JAC <- .lavJacobianD(func=fit@Model@def.function, x=x)
  ## Compute the sampling covariance matrix using delta method
  V <- JAC %*% lavaan::inspect(fit, what="vcov") %*% t(JAC)
  ## Add the variable names for ease of reference
  dimnames(V) <- list(names(y), names(y))
  out <- list(y=y, V=V)
  }
  
  out
}
