osmasem3L <- function(model.name="osmasem3L", RAM=NULL, cluster=NULL,
                      RE.typeB=c("Diag", "Symm"), RE.typeW=c("Diag", "Symm"),
                      data, Mmatrix=NULL, TmatrixB=NULL, TmatrixW=NULL,
                      Ax=NULL, Sx=NULL, A.lbound=NULL, A.ubound=NULL,
                      subset.variables=NULL, subset.rows=NULL, 
                      intervals.type = c("z", "LB"),
                      mxModel.Args=NULL, mxRun.Args=NULL,
                      suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {

    RE.typeB <- match.arg(RE.typeB)
    RE.typeW <- match.arg(RE.typeW)
    
    intervals.type <- match.arg(intervals.type)
    switch(intervals.type,
           z = intervals <- FALSE,
           LB = intervals <- TRUE)

    ## Operator for not in
    ## "%ni%" <- Negate("%in%")

    ## Filter out variables not used in the analysis
    if (!is.null(subset.variables)) {
        index  <- !(subset.variables %in% data$obslabels)
        if (any(index)) {
            stop(paste0(paste0(subset.variables[index], collapse = ", "), " are not in the ylabels.\n"))
        }
        temp <- .genCorNames(subset.variables)
        obslabels <- subset.variables
        ylabels <- temp$ylabels
        vlabels <- temp$vlabels        
    } else {
        obslabels <- data$obslabels
        ylabels <- data$ylabels
        vlabels <- data$vlabels
    }

    ## no. of y (effect sizes)
    p <- length(ylabels)

    ## If RAM is provided, create the matrices based on it.
    if (!is.null(RAM)) {
        Mmatrix <- create.vechsR(A0=RAM$A, S0=RAM$S, F0=RAM$F, Ax=Ax, Sx=Sx,
                                 A.lbound=A.lbound, A.ubound=A.ubound)
        TmatrixB <- create.Tau2(RAM=RAM, RE.type=RE.typeB, Transform="expLog", level="between")
        TmatrixW <- create.Tau2(RAM=RAM, RE.type=RE.typeW, Transform="expLog", level="within")
    }

    ## Get the labels of the latent variables in the variance component (between)
    latlabelsB <- TmatrixB$vecTau1$labels
    
    ## Create known sampling variance covariance matrix
    Vmatrix <- create.V(vlabels, type="Symm", as.mxMatrix=TRUE)
   
    if (is.null(subset.rows)) {
        subset.rows <- rep(TRUE, nrow(data$data))
    }
    
    ## Select the within data
    mx.data <- data$data[subset.rows, ]
    mx.data[, cluster] <- as.factor(mx.data[, cluster])

    ## Create the between data
    mx.dataB <- eval(parse(text=paste0("data.frame(", cluster,
                                       "=as.factor(unique(mx.data[, '", cluster, "'])))")))

    ## modelB <- mxModel(model='Bet', type='RAM', latentVars=latlabelsB, TmatrixB,
    ##                        mxData(observed=mx.dataB, type='raw', primaryKey=cluster),
    ##                        mxMatrix('Zero', p, p, name='A', dimnames=list(latlabelsB, latlabelsB)),
    ##                        mxMatrix('Zero', 0, p, name='F', dimnames=list(NULL, latlabelsB)),
    ##                        mxMatrix('Zero', 1, p, name='M', dimnames=list(NULL, latlabelsB)),
    ##                        mxExpectationRAM('A', 'Tau2B', 'F', 'M'))

    ## Between model
    ## A: a pxp zero matrix
    ## S: a pxp Tau2 (between)
    ## M: a 1xp zero vector
    ## Trick to avoid the warning of "no visible binding for global variable"
    modelB <- eval(parse(text="mxModel(model='B', type='RAM', latentVars=latlabelsB, TmatrixB,
                           mxData(observed=mx.dataB, type='raw', primaryKey=cluster),
                           mxMatrix('Zero', p, p, name='A', dimnames=list(latlabelsB, latlabelsB)),
                           mxMatrix('Zero', 0, p, name='F', dimnames=list(NULL, latlabelsB)),
                           mxMatrix('Zero', 1, p, name='M', dimnames=list(NULL, latlabelsB)),
                           mxExpectationRAM(A='A', S='Tau2B', F='F', M='M'))"))

    ## Need to add the dimnames in the model-implied mean structure in the Within model
    Mmatrix$vechsR$.dimnames <- list(NULL, ylabels)
    
    ## mx.model <- mxModel(model=model.name, type='RAM', modelB, Mmatrix, TmatrixW, Vmatrix,
    ##                     manifestVars=ylabels,
    ##                     mxData(mx.data, 'raw'),
    ##                     mxAlgebra(Tau2W+V, name='expCov'),
    ##                     mxMatrix('Zero', p, p, name='A', dimnames=list(ylabels, ylabels)),
    ##                     mxMatrix('Iden', p, p, name='F', dimnames=list(ylabels, ylabels)),
    ##                     mxMatrix('Full', p, p, FALSE, 1, name = 'T',
    ##                              joinKey=cluster, joinModel='Bet', dimnames=list(ylabels, latlabelsB)),
    ##                     mxExpectationRAM('A', 'expCov', 'F', 'vechsR', between='T'),
    ##                     mxCI(c('Amatrix', 'Smatrix', 'Tau2W', 'Bet.Tau2B')))

    ## Within model
    ## A: a pxp zero matrix
    ## S: a pxp Tau2 (within + V)
    ## M: a 1xp vector of model-implied R    
    mx.model <- eval(parse(text="mxModel(model=model.name, type='RAM', modelB,
                         Mmatrix, TmatrixW, Vmatrix, manifestVars=ylabels,
                         mxData(mx.data, 'raw'),
                         mxAlgebra(Tau2W+V, name='expCov'),
                         mxMatrix('Zero', p, p, name='A', dimnames=list(ylabels, ylabels)),
                         mxMatrix('Iden', p, p, name='F', dimnames=list(ylabels, ylabels)),
                         mxMatrix('Full', p, p, FALSE, 1, name = 'T',
                                  joinKey=cluster, joinModel='B', dimnames=list(ylabels, latlabelsB)),
                         mxExpectationRAM(A='A', S='expCov', F='F', M='vechsR', between='T'),
                         mxCI(c('Amatrix', 'Smatrix', 'Tau2W', 'Bet.Tau2B')))"))
   
    ## Add additiona arguments from RAM to mxModel
    if (!is.null(RAM$mxalgebras)) {
        for (i in seq_along(RAM$mxalgebras)) {
            mx.model <- mxModel(mx.model, RAM$mxalgebras[[i]])
        }
        ## check if they are mxalgebra, not mxconstraint
        algebra.names <- names(RAM$mxalgebras)
        isalgebra <- !grepl("^constraint[0-9]+", algebra.names)
        if (any(isalgebra)) {
            mx.model <- mxModel(mx.model, mxCI(algebra.names[isalgebra]))
        }
    }  

    ## Add additional arguments to mxModel
    if (!is.null(mxModel.Args)) {
        for (i in seq_along(mxModel.Args)) {
            mx.model <- mxModel(mx.model, mxModel.Args[[i]])
        }
    }

    ## Assign same starting values, e.g., A and S matrices are different
    mx.model <- omxAssignFirstParameters(mx.model)
    
    ## Return mx model without running the analysis
    if (run==FALSE) {
        return(mx.model)
    } else {        
        mx.fit <- tryCatch(do.call(mxRun, c(list(mx.model, intervals=intervals,
                                                 suppressWarnings=suppressWarnings,
                                                 silent=silent),
                                            mxRun.Args)), error = function(e) e)
    }

    # try to run it with error message as output
    if (inherits(mx.fit, "error")) {
        warning(print(mx.fit))
    }

    out <- list(call=match.call(), 
                Mmatrix=Mmatrix, 
                TmatrixB=TmatrixB,
                TmatrixW=TmatrixW,
                Vmatrix=Vmatrix, 
                data=data,
                cluster=cluster,
                labels=list(obslabels=obslabels, ylabels=ylabels, vlabels=vlabels),
                mxModel.Args=mxModel.Args,
                mxRun.Args=mxRun.Args,
                subset.variables=subset.variables,
                subset.rows=subset.rows,
                mx.model=mx.model, mx.fit=mx.fit)
    class(out) <- 'osmasem3L'
    out
}

coef.osmasem3L <- function(object, select=c("fixed", "all", "random"), ...) {
    if (!is.element("osmasem3L", class(object)))
        stop("\"object\" must be an object of class \"osmasem3L\".")
    
    mx.coef <- coef(object$mx.fit)
    ## index for untransformed random effects (not the correct ones!) 
    index <- grep("Tau1_|Cor_", names(mx.coef))

    select <- match.arg(select)
    switch(select,
           ## all = my.name <- my.name,
           fixed =  mx.coef <- mx.coef[-index],
           random = mx.coef <- mx.coef[index])

    if (select!="fixed") {
        warning("\"Tau1_xx\" is not the variance component of the random effects.\nPlease use VarCorr() to get the variance component.\n")
    }
    
    ## Rearrange the elements according to the original variable order if it is a tssem3L1 object.
    if (is.element("tssem3L1", class(object))) {
        mx.coef[object$data$ylabels]
    } else {    
        mx.coef
    }
}

vcov.osmasem3L <- function(object, select=c("fixed", "all", "random"), ...) {
    if (!is.element("osmasem3L", class(object)))
        stop("\"object\" must be an object of class \"osmasem3L\".")

    # labels of the parameters    
    ## my.name <- summary(object$mx.fit)$parameters$name
    my.name <- names(omxGetParameters(object$mx.fit))
    my.name <- my.name[!is.na(my.name)]

    ## index for untransformed random effects (not the correct ones!) 
    index <- grep("Tau1_|Cor_", my.name) 
    
    select <- match.arg(select)
    switch( select,
         ## all = my.name <- my.name,
         fixed =  my.name <- my.name[-index],
         random = my.name <- my.name[index])
    
    if (select!="fixed") {
        warning("\"Tau1_xx\" is not the variance component of the random effects.")
    }

    out <- vcov(object$mx.fit)
    out <- out[my.name, my.name, drop=FALSE]

    ## Rearrange the elements according to the original variable order if it is a tssem3L1 object.
    if (is.element("tssem3L1", class(object))) {
        out[object$data$ylabels, object$data$ylabels, drop=FALSE]
    } else {    
        out
    }       
}


## Fit a saturated model with either a diagonal or symmetric variance component of random effects
.osmasemSatIndMod3L <- function(osmasem.obj=NULL, model=c("Saturated", "Independence"),
                                Std.Error=FALSE, extraTries=50, run=TRUE, ...) {

    if (!is.element("osmasem3L", class(osmasem.obj)))
        stop("\"osmasem.obj\" must be an object of class \"osmasem3L\".")

    TmatrixB <- osmasem.obj$TmatrixB
    TmatrixW <- osmasem.obj$TmatrixW
    Vmatrix <- osmasem.obj$Vmatrix

    ## No. of variables in the model
    p <- nrow(TmatrixB$Cor$values)
    data <- osmasem.obj$data

    model <- match.arg(model)
    cluster <- osmasem.obj$cluster
    ylabels <- osmasem.obj$labels$ylabels

    if (model=="Saturated") {
        ## Saturated model
        Mmatrix <- mxMatrix(type="Full", free=TRUE, labels=paste0("Mu", seq_len(p)),
                            nrow=1, ncol=p, dimnames=list(NULL, ylabels), name="Mu")
    } else {
        ## Independence model
        Mmatrix <- mxMatrix(type="Full", free=FALSE, labels=paste0("Mu", seq_len(p)),
                            nrow=1, ncol=p, , dimnames=list(NULL, ylabels), name="Mu")        
    }
    
    ## Select data for the analaysis
    mx.data <- data$data[osmasem.obj$subset.rows, ]
    mx.data[, cluster] <- as.factor(mx.data[, cluster])
    
    ## Get the labels of the latent variables in the variance component (between)
    latlabelsB <- TmatrixB$vecTau1$labels

    ## Between data
    mx.dataB <- eval(parse(text=paste0("data.frame(", cluster,
                                       "=as.factor(unique(mx.data[, '", cluster, "'])))")))

    ## Trick to avoid the warning of "no visible binding for global variable"
    modelB <- eval(parse(text="mxModel(model='B', type='RAM', latentVars=latlabelsB, TmatrixB,
                           mxData(observed=mx.dataB, type='raw', primaryKey=cluster),
                           mxMatrix('Zero', p, p, name='A', dimnames=list(latlabelsB, latlabelsB)),
                           mxMatrix('Zero', 0, p, name='F', dimnames=list(NULL, latlabelsB)),
                           mxMatrix('Zero', 1, p, name='M', dimnames=list(NULL, latlabelsB)),
                           mxExpectationRAM('A', 'Tau2B', 'F', 'M'))"))
    
    mx.model <- eval(parse(text="mxModel(model=model, type='RAM', modelB,
                          Mmatrix, TmatrixW, Vmatrix, manifestVars=ylabels,
                          mxData(mx.data, 'raw'),
                          mxAlgebra(Tau2W+V, name='expCov'),
                          mxMatrix('Zero', p, p, name='A', dimnames=list(ylabels, ylabels)),
                          mxMatrix('Iden', p, p, name='F', dimnames=list(ylabels, ylabels)),
                          mxMatrix('Full', p, p, FALSE, 1, name = 'T',
                                   joinKey=cluster, joinModel='B', dimnames=list(ylabels, latlabelsB)),
                          mxExpectationRAM('A', 'expCov', 'F', 'Mu', between='T'))"))
     
    ## It may speed up the analysis
    if (Std.Error==FALSE) {
        mx.model <- mxOption(mx.model, "Calculate Hessian", "No")
        mx.model <- mxOption(mx.model, "Standard Errors", "No")
    }
    
    ## Return mx model without running the analysis
    if (run==FALSE) {
        return(mx.model)
    } else {
        mx.fit <- mxTryHard(mx.model, extraTries=extraTries, silent=TRUE, ...) 
    }
    
    # try to run it with error message as output
    ## if (inherits(mx.fit, "error")) {
    ##     warning(print(mx.fit))
    ## }
    mx.fit
}

anova.osmasem3L <- function(object, ..., all=FALSE) {
  base <- lapply(list(object), function(x) x$mx.fit)
  comparison <- lapply(list(...), function(x) x$mx.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

summary.osmasem3L <- function(object, fitIndices=FALSE, numObs, ...) {
    if (!is.element("osmasem3L", class(object)))
        stop("\"object\" must be an object of class \"osmasem3L\".")

    ## If numObs is not provided, use the total N
    if (missing(numObs)) {
        ## numObs <- sum(object$data$n)-length(object$data$n)
        numObs <- sum(object$data$n[object$subset.rows])
    }

    ## Calculate chi-square statistic and other fit indices
    if (fitIndices) {
        Sat.stat <- .osmasemSatIndMod3L(object, model="Saturated", Std.Error=FALSE)
        Ind.stat <- .osmasemSatIndMod3L(object, model="Independence", Std.Error=FALSE) 
    
        ## Sat.model=FALSE for saturated model if there are either errors or nonconvergent.
        if ( inherits(Sat.stat, "error") | !(Sat.stat$output$status$code %in% c(0,1)) ) {
            Sat.model <- FALSE
        } else {
            Sat.stat <- summary(Sat.stat)
            if (Sat.stat$degreesOfFreedom <= 0) {
                Sat.model <- FALSE
            } else {
                Sat.model <- TRUE
            }
        }

        ## Ind.model=FALSE for independence model if there are either errors or nonconvergent.
        ## NA is acceptable for independence model as it means that optimization was not attempted.
        if ( inherits(Ind.stat, "error") | !(Ind.stat$output$status$code %in% c(NA, 0, 1)) ) {
            Ind.model <- FALSE
        } else {
            Ind.stat <- summary(Ind.stat)
            if (Ind.stat$degreesOfFreedom <= 0) {
                Ind.model <- FALSE
            } else {
                Ind.model <- TRUE
            }
        }

        ## If sat.model is defined
        if (Sat.model) {
            if (Ind.model) {
                out <- summary(object$mx.fit,
                               SaturatedLikelihood=Sat.stat$Minus2LogLikelihood,
                               SaturatedDoF=Sat.stat$degreesOfFreedom,
                               IndependenceLikelihood=Ind.stat$Minus2LogLikelihood,
                               IndependenceDoF=Ind.stat$degreesOfFreedom,
                               numObs=numObs, ...)
            } else {
                out <- summary(object$mx.fit,
                               SaturatedLikelihood=Sat.stat$Minus2LogLikelihood,
                               SaturatedDoF=Sat.stat$degreesOfFreedom,
                               numObs=numObs, ...)
                warning("There are errors in either fitting the independence model or its degree of freedom being non-positive.\n")
        }
            ## if Sat.model is not defined, no chi-square and fit indices    
        } else {
            out <- summary(object$mx.fit, numObs=numObs, ...)
            warning("There are errors in either fitting the saturated model or its degree of freedom being non-positive.\n")
        }

        ## fitIndices=FALSE
    } else {
        out <- summary(object$mx.fit, numObs=numObs, ...)
    }

    ## Additional output to the mx summary
    out$parameters$`z value` <- with(out$parameters, Estimate/Std.Error)
    out$parameters$`Pr(>|z|)` <- 2*(1-pnorm(abs(out$parameters$`z value`)))
    out
}

osmasem3LR2 <- function(model1, model0, R2.truncate=TRUE) {
    if (!all(c(class(model0), class(model1)) %in% "osmasem3L"))
        stop("Both \"model0\" and \"model1\" must be objects of class \"osmasem3L\".")

    ## Model 0
    Tau2.0 <- VarCorr(model0)
    ## Within
    Tau2.0W <- diag(Tau2.0[["Tau2W"]])
    ## Between
    Tau2.0B <- diag(Tau2.0[["Tau2B"]])

    ## Model 1
    Tau2.1 <- VarCorr(model1)
    ## Within
    Tau2.1W <- diag(Tau2.1[["Tau2W"]])
    ## Between
    Tau2.1B <- diag(Tau2.1[["Tau2B"]])   
    
    ## Within
    R2W <- (Tau2.0W-Tau2.1W)/Tau2.0W

    ## Between
    R2B <- (Tau2.0B-Tau2.1B)/Tau2.0B    

    ## Truncate negative R2 to 0
    if (R2.truncate) {
        R2W[R2W<0] <- 0
        R2B[R2B<0] <- 0
    }
    
    p <- nrow(model0$TmatrixW$Cor$values)        
    names.tau2 <- paste0("Tau2_", seq_len(p), "_", seq_len(p))
    names(R2W) <- names(Tau2.0W) <- names(Tau2.1W) <- paste0(names.tau2, "W")
    names(R2B) <- names(Tau2.0B) <- names(Tau2.1B) <- paste0(names.tau2, "B")
    
    list(Tau2.0W=Tau2.0W, Tau2.1W=Tau2.1W, R2W=R2W,
         Tau2.0B=Tau2.0B, Tau2.1B=Tau2.1B, R2B=R2B)
}

osmasem3LSRMR <- function(x) {
    ## Check if there are moderators in either A1 or S1
    if (!is.null(x$Mmatrix$A1) | !is.null(x$Mmatrix$S1))
        stop("Moderators are not allowed in calculating the SRMR in osmasem3L.\n")

    ## Saturated model which should be very close to the stage 1 results in the TSSEM
    Sat.stat <- .osmasemSatIndMod3L(x, model="Saturated", Std.Error=FALSE)

    ## Similar to the sample correlation matrix
    sampleR <- vec2symMat(eval(parse(text = "mxEval(Mu, Sat.stat)")), diag=FALSE)

    ## Model implied correlation matrix
    impliedR <- mxEval(impliedR, x$mx.fit)

    ## SRMR (Cheung, 2015 book, Eq. 2.36)
    sqrt(mean(vechs(sampleR-impliedR)^2))
}


