#' Convert \code{lavaan} models to RAM models
#' 
#' It converts models specified in \code{lavaan} model syntax to RAM models.
#' 
#' It uses the \code{\link[lavaan]{model.syntax}} to do the conversion.
#' Experimental: functions of parameters (`:=` in lavaan) and constraints
#' (`==`, `>`, and `<` in lavaan) will be converted to mxAlgebra and
#' mxConstraint in OpenMx. As there are differences between lavaan and OpenMx,
#' they may not work properly.
#' 
#' @param model A character string of model using the lavaan model syntax.
#' @param obs.variables A character vector of the observed variables. The
#' observed variables in the RAM specification will follow the order specified
#' in \code{obs.variables}. It is important to check whether the order of the
#' observed variables matches the order in the dataset.
#' @param A.notation A character string to be used in the A matrix if the
#' labels are not included in the lavaan model. For example, the label will be
#' "yONx" for regressing "y" on "x".
#' @param S.notation A character string to be used in the S matrix if the
#' labels are not included in the lavaan model. For example, the label will be
#' "yWITHx" for the covariance between "y" with "x" and "yWITHy" for the
#' (error) variance of "y".
#' @param M.notation A character string to be used in the M matrix if the
#' labels are not included in the lavaan model. For example, the label will be
#' "ymean" for the mean of "y" if \code{M.notation="mean"}.
#' @param A.start A numeric value of starting value for the \code{Amatrix} when
#' the starting values are not provided.
#' @param S.start A numeric value of starting value for the \code{Smatrix} when
#' the starting values are not provided.
#' @param M.start A numeric value of starting value for the \code{Mmatrix} when
#' the starting values are not provided.
#' @param auto.var Logical. If \code{TRUE}, the residual variances and the
#' variances of exogenous latent variables are included in the model and set
#' free. See \code{\link[lavaan]{model.syntax}}.
#' @param std.lv Logical. If \code{TRUE}, the metric of each latent variable is
#' determined by fixing their variances to 1.0. If FALSE, the metric of each
#' latent variable is determined by fixing the factor loading of the first
#' indicator to 1.0. See \code{\link[lavaan]{model.syntax}}.
#' @param ngroups Number. The number of groups in the \code{model}. See
#' \code{\link[lavaan]{model.syntax}}.
#' @param \dots Further arguments to be passed to
#' \code{\link[lavaan]{model.syntax}}
#' @return A list of RAM specification with \code{A}, \code{S}, \code{F}, and
#' \code{M} matrices.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[semPlot]{ramModel}}, \code{\link[metaSEM]{Becker92}},
#' \code{\link[metaSEM]{Becker09}}, \code{\link[metaSEM]{Digman97}},
#' \code{\link[metaSEM]{Hunter83}}, \code{\link[metaSEM]{as.mxMatrix}},
#' \code{\link[metaSEM]{checkRAM}}
#' @keywords methods tssem
#' @examples
#' 
#' ## Regression model on correlation matrix
#' model1 <- "## y is modelled by x1, x2, and x3
#'            y ~ b1*x1 + b2*x2 + b3*x3
#'            ## Fix the independent variables at 1
#'            x1 ~~ 1*x1
#'            x2 ~~ 1*x2
#'            x3 ~~ 1*x3
#'            ## Declare the correlations among the independent variables
#'            x1 ~~ x2
#'            x1 ~~ x3
#'            x2 ~~ x3
#'            ## Constraint
#'            b3 == b1 + b2
#'            ## Function of parameters
#'            fn1 := b1*b2^b3"
#' 
#' ## Compare the arrangements of variables with and without
#' ## specifying the obs.variables arguments. 
#' lavaan2RAM(model1, obs.variables=c("y", "x1", "x2", "x3"))
#' 
#' ## Two-factor CFA model
#' model2 <- "f1 =~ x1 + x2 + x3
#'            f2 =~ x4 + x5 + x6
#'            ## Declare the correlation between f1 and f2
#'            ## and label it with cor_f1f2
#'            f1 ~~ cor_f1f2*f2"
#' 
#' lavaan2RAM(model2)
#' 
#' ## Regression model with the mean structure
#' model3 <- "y ~ x
#'            ## Intercept of y
#'            y ~ 1
#'            ## Mean of x
#'            x ~ 1"
#' 
#' lavaan2RAM(model3)
#' 
#' ## Multiple group regression model
#' ## Different intercepts with a common slope
#' model4 <- "y ~ c(a1, a2)*1 + c(b, b)*x"
#' 
#' lavaan2RAM(model4, ngroups=2)
#' 
lavaan2RAM <- function(model, obs.variables = NULL, A.notation="ON", S.notation="WITH",
                       M.notation="mean", A.start=0.1, S.start=0.5, M.start=0,
                       auto.var = TRUE, std.lv = TRUE,
                       ngroups=1, ...) {
  ## if (!requireNamespace("lavaan", quietly=TRUE))    
  ##     stop("\"lavaan\" package is required for this function.")

  ## Default: fix the latent independent variables at 1
  my.model <- lavaan::lavaanify(model, fixed.x = FALSE, auto.var=auto.var,
                                std.lv=std.lv, ngroups=ngroups, ...)

  ## Maximum no. of groups
  max.gp <- max(my.model$group)

  ## Empty list to store the matrices per group
  out <- list()
    
  for (gp in seq_len(max.gp)) {
    ## Select the ith group
    mod <- my.model[my.model$group==gp, ]

    ## if (any...) is required to avoid error when there is no element for the assignment
    ## set the starting values in A if NA
    if (any((mod$op=="=~"|mod$op=="~")&is.na(mod$ustart))) {            
      mod[(mod$op=="=~"|mod$op=="~")&is.na(mod$ustart), ]$ustart <- A.start
    }
        
    ## set the starting values in M if NA
    if (any(mod$op=="~1"&is.na(mod$ustart))) {
      mod[mod$op=="~1"&is.na(mod$ustart), ]$ustart <- M.start      
    }
        
    ## set the starting values in S and free parameters if NA (variances)
    if (any(mod$op=="~~"&is.na(mod$ustart)&(mod$lhs==mod$rhs))) {
      mod[mod$op=="~~"&is.na(mod$ustart)&(mod$lhs==mod$rhs), ]$ustart <- S.start
    }
        
    ## Set the starting values in S and free parameters if NA (covariances)
    if (any(mod$op=="~~"&is.na(mod$ustart)&(mod$lhs!=mod$rhs))) {
      mod[mod$op=="~~"&is.na(mod$ustart)&(mod$lhs!=mod$rhs), ]$ustart <- 0
    }

    ## all variables
    ## Removed sort(); otherwise, the variables will be arranged as x1, x10, x2, x3...
    all.var <- unique(c(mod$lhs, mod$rhs))
        
    ## latent variables: (with indicators)
    latent <- unique(mod[mod$op== "=~", ]$lhs)
    ## observed variables: not latent
    observed <- all.var[!(all.var %in% latent)]
    ## remove empty string "" when there are mean structure
    observed <- observed[observed !=""]
        
    ## check whether observed in model = observed in argument
    if (!is.null(obs.variables)) {
      if (!identical(sort(observed), sort(obs.variables))) {
        stop("Names in \"obs.variables\" do not agree with those in model.\n")
      } else {
        ## arrange the observed variables according to obs.var argument
        observed <- obs.variables
      }
    }

    ## if there are latent variables
    if (length(latent)>0) {
      ## arrange variable list as observed + latent
      all.var <- c(observed, latent)
    } else {
      all.var <- observed
    }
    
    no.lat <- length(latent)
    no.obs <- length(observed)
    no.all <- no.lat+no.obs
    
    Amatrix <- matrix(0, ncol=no.all, nrow=no.all, dimnames=list(all.var, all.var))
    Smatrix <- matrix(0, ncol=no.all, nrow=no.all, dimnames=list(all.var, all.var))
    ## Fixed Mmatrix by setting the default to estimate the means;
    ## otherwise, the mean structure is misspecified.
    ## Defaults: observed variables are free, whereas latent variables are fixed at 0
    Mmatrix <- matrix(c(paste0("0*", observed, "mean"), rep("0", no.lat)),
                      nrow=1, ncol=no.all, dimnames=list(1, all.var))

    ## Prepare the labels
    for (i in seq_len(nrow(mod))) {
      ## if there is no label
      if (mod[i, ]$label=="") {
        switch(mod[i, ]$op,
               "=~" = mod[i, ]$label <- paste0(mod[i, ]$rhs, A.notation, mod[i, ]$lhs),
               "<~" = mod[i, ]$label <- paste0(mod[i, ]$lhs, A.notation, mod[i, ]$rhs),
               "~"  = mod[i, ]$label <- paste0(mod[i, ]$lhs, A.notation, mod[i, ]$rhs),
               "~~" = mod[i, ]$label <- paste0(mod[i, ]$lhs, S.notation, mod[i, ]$rhs),
               "~1" = mod[i, ]$label <- paste0(mod[i, ]$lhs, M.notation))
            }
    }
    
    ## Replace NA to 0 in ustart if there are still NA
    ## mod$ustart[is.na(mod$ustart)] <- 0
        
    ## keys in as.mxMatrix format
    key <- with(mod, ifelse(free==0, yes=ustart, no=paste(ustart, label, sep="*")))  
        
    for (i in seq_len(nrow(mod))) {
      my.line <- mod[i, ]
      switch(my.line$op,
             ## lhs: IV; rhs: DV
             "=~" = Amatrix[my.line$rhs, my.line$lhs] <- key[i],
             ## lhs: DV; rhs: IV
             "<~" = Amatrix[my.line$lhs, my.line$rhs] <- key[i],
             "~"  = Amatrix[my.line$lhs, my.line$rhs] <- key[i],
             "~~" = Smatrix[my.line$lhs, my.line$rhs] <-
               Smatrix[my.line$rhs, my.line$lhs] <- key[i],
             ## means
             "~1" = Mmatrix[1, my.line$lhs] <- key[i]
             )  ## from switch
    }  ## from for loop

    Fmatrix <- create.Fmatrix(c(rep(1, no.obs), rep(0, no.lat)), as.mxMatrix=FALSE)
    dimnames(Fmatrix) <- list(observed, all.var)
    
    out[[gp]] <- list(A=Amatrix, S=Smatrix, F=Fmatrix, M=Mmatrix)
  }

  ## Add group names, 1, 2, 3... to the list
  names(out) <- seq_along(out)
    
  ## If there are constraints such as .p1.==.p2.; remove them first
  ## otherwise, .p1.==.p2. will create an empty list in mxalgebra
  if (length(grep("^\\.", my.model$lhs)) >0 ) {
    my.model <- my.model[-grep("^\\.", my.model$lhs), ]
  }
    
  ## Check if there are constraints or algebras
  if (any(my.model$group==0)) {
        
    ## An empty list to store mxAlgebra and mxConstraint
    mxalgebra <- list()
    con_index <- 1   # A counter for no. of constraints, started from 1

    ## Constraints and algebras only
    y <- my.model[my.model$group==0, , drop=FALSE]
    
    for (i in seq_len(nrow(y))) {
      switch(y[i, 'op'],
             ## mxAlgebra
             ":=" = { eval(parse(text=paste0(y[i,'lhs'],
                                             "<- mxAlgebra(", y[i,'rhs'],
                                             ", name=\"", y[i,'lhs'], "\")")))
                                             eval(parse(text=paste0("mxalgebra <-
                                                  c(mxalgebra, ", y[i, 'lhs'],
                                                  "=", y[i, 'lhs'], ")")))  },
             ## Default condition to test if there are mxConstraints
             if (y[i, 'op'] %in% c("==", ">", "<")) {
               eval(parse(text=paste0("constraint", con_index, " <- mxConstraint(", y[i, 'lhs'],
                                      y[i, 'op'], y[i, 'rhs'], ", name=\"constraint", con_index, "\")")))
               eval(parse(text=paste0("mxalgebra <- c(mxalgebra, constraint",
                                      con_index, "=constraint", con_index, ")")))
               con_index <- con_index + 1
             }
             )  ## from switch
    }  ## from for loop

    ## Append mxalgebra to out[[1]]
    out[[1]]  <- list(A=out[[1]]$A, S=out[[1]]$S, F=out[[1]]$F, M=out[[1]]$M,
                      mxalgebras=mxalgebra)
  }  ## End if
  ## else {
  ##     out[[1]]  <- list(A=out[[1]]$A, S=out[[1]]$S, F=out[[1]]$F, M=out[[1]]$M)
  ## }


  ## Output the first list instead of a list of one item when there is only 1 group
  if (max.gp==1) {
    out  <- out[[1]]
  }
    
  out
}
