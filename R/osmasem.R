Cor2DataFrame <- function(x, n, v.na.replace=TRUE, row.names.unique=FALSE,
                          cor.analysis=TRUE, ...) {
    if (length(x) != length(n)) stop("Lengths of 'x' and 'n' are different.\n")
    
    if (cor.analysis) {
        my.df <- list2matrix(x=suppressWarnings(lapply(x, cov2cor)), diag=FALSE)
    } else {
        my.df <- list2matrix(x=x, diag=TRUE)
    }

    acovR <- asyCov(x=x, n=n, cor.analysis=cor.analysis, ...)

    ## NA is not allowed in definition variables
    ## They are replaced by 1e10
    if (v.na.replace) acovR[is.na(acovR)] <- 1e10

    data <- suppressWarnings(data.frame(my.df, acovR, check.names=FALSE))
    
    ## Use unique row names if the row names are duplicated.
    if (row.names.unique) rownames(data) <- make.names(names(x), unique=TRUE)    

    list(data=data, n=n, ylabels=dimnames(my.df)[[2]], vlabels=dimnames(acovR)[[2]])
}


## A function to extract levels of path directions started from IVs to DVs
.pathLevels <- function(A) {
    ## assuming non-recursive models
    A <- A$free
        
    ## no. of variables
    p <- ncol(A)
    dv <- apply(A, 1, any)
    ## iv: variables pointed towards other variables
    iv <- apply(A, 2, any)

    Level <- list()
    ## IVs
    Level[[1]] <- iv==TRUE&dv==FALSE
    ## no. of variables counted
    count <- length((1:p)[Level[[1]]])

    ## do until all p variables are countered
    while(p > count) {
        level <- length(Level)

        ## predicated by previous level
        cand1 <- A[, Level[[level]], drop=FALSE]
        cand1 <- apply(cand1, 1, any)

        ## predicted by cand1
        cand2 <- A[, cand1, drop=FALSE]
        cand2 <- apply(cand2, 1, any)
        Level[[level+1]] <- cand1==TRUE&cand2==FALSE

        ## no. of elements counted
        count <- count + length((1:p)[Level[[level+1]]])
    }

    list(IV=Level[1], DV=Level[-1])
}



## A0: Intercept
## A1: 1st moderator
## Amatrix = A0 + A1*Ax[[1]] + A1*Ax[[2]]
## Smatrix = S0 + S1*Sx[[1]] + S2*Sx[[2]]
create.vechsR <- function(A0, S0, F0=NULL, Ax=NULL, Sx=NULL) {

    if (is.matrix(A0)) {
        A0 <- as.mxMatrix(A0, name="A0")
    } else {
        A0@name <- "A0"
    }    

    if (is.matrix(S0)) {
        S0 <- as.mxMatrix(S0, name="S0")
    } else {
        S0@name <- "S0"
    }

    ## Basic checking
    checkRAM(Amatrix=A0, Smatrix=S0, cor.analysis=TRUE)
    ## ## Additional check on whether the diagonals of S0 are 0
    ## if (any(Diag(S0$free)==TRUE)) {
    ##     stop("Diagonal of 'S0' must be 0.\n")
    ## }       
    
    ## Create an Identity matrix
    Iden <- create.Fmatrix(rep(1, ncol(A0$values)), name="Iden")

    ## If F0 is not specified, assume all variables are selected.
    if (is.null(F0)) {
        Fmatrix <- Iden
        Fmatrix$name <- "Fmatrix"
    } else {
        Fmatrix <- as.mxMatrix(F0, name="Fmatrix")
    }

    ## Consider error variances as functions of parameters
    ## Check types of variables
    ## dv (not ivs): variables pointed by some variables
    dv <- apply(A0$free, 1, any)

    ## ## iv: variables pointed towards other variables
    ## iv <- apply(Amatrix$free, 2, any)
    ## ## med: mediators
    ## med <- iv & dv

    ## Error variances of dv are modified as 0.
    ## S0: Modified Smatrix without error variances on the dependent variables
    ## Fix the diagonal elements associated with dv to 0
    diag(S0$labels)[dv] <- NA
    diag(S0$values)[dv] <- 0
    diag(S0$free)[dv] <- FALSE

    ## ## Use it rather than Smatrix in mxCI()
    ## S1_labels <- S1$labels
    ## diag(S1_labels) <- NA
    ## S1_labels <- c(na.omit(unique(c(S1_labels))))

    ## Levels of directions, which variables are dv
    Level <- .pathLevels(A0)$DV

    ## When there is no predictor, only A0.
    ## Amatrix=A0
    ## ToCheck: Use mxalgebra or copy A0 to Amatrix?
    if (is.null(Ax)) {
        Amatrix <- mxAlgebra(A0, name="Amatrix")
    } else {
        ## Convert Apred into a list if it is not
        if (!is.list(Ax)) Ax <- list(Ax)

        ## Basic checking on Ax
        lapply(Ax, function(x) checkRAM(Amatrix=x))
        
        for (i in seq_len(length(Ax))) {

            ## Parameters for the moderator based on the Ax, not A0
            ## When the labels include "data", they are labels of the moderator.
            index_mod <- apply(Ax[[i]], 1:2, function(x) grepl("data", x))

            ## A temp A0. The free parameters are based on Ax, not A0
            Atemp <- A0
            Atemp$free <- index_mod
            Atemp$labels[!index_mod] <- NA

            text1 <- paste0("A",i, " <- Atemp; A",i, "$name <- 'A",i, "'")
            eval(parse(text=text1))

            ## Labels = Labels in A0 + "_1", "_2", etc.
            text2 <- paste0("A",i, "$labels <- apply(A",i,
                            "$labels, 1:2, function(x) {if (is.na(x)) NA else paste0(x,", "'_",i, "')})")
            eval(parse(text=text2))            

            ## Values of the moderators with definition variables
            text3 <- paste0("Ax",i, " <- as.mxMatrix(Ax[[",
                            i, "]], name='Ax",i, "')")
            eval(parse(text=text3))
        }  ## from 'for'

        ## Amatrix= A0+A1*Ax1+A2*Ax2+...
        text4 <- paste0("A", seq_len(length(Ax)), "*Ax",  seq_len(length(Ax)),
                        collapse=" + ")   
        text5 <- paste0("Amatrix <- mxAlgebra(A0 + ", text4, ", name='Amatrix')")
        eval(parse(text=text5))
        
    }   ## from 'if else' 

    ## When there is no predictor, only S0.
    ## SnoErr=S0
    if (is.null(Sx)) {
        SnoErr <- mxAlgebra(S0, name="SnoErr")
    } else {
        ## Convert Sx into a list if it is not
        if (!is.list(Sx)) Sx <- list(Sx)

        ## Basic checking on Sx
        ## Don't check the diagonals
        lapply(Sx, function(x) checkRAM(Smatrix=x, cor.analysis=FALSE))
        
        for (i in seq_len(length(Sx))) {

            ## Parameters for the moderator
            text1 <- paste0("S",i, " <- S0; S",i, "$name <- 'S",i, "'")
            eval(parse(text=text1))

            ## Relabel the parameter labels
            ## Labels = Labels in S0 + "_1", "_2", etc.
            text2 <- paste0("S",i, "$labels <- apply(S",i,
                            "$labels, 1:2, function(x) {if (is.na(x)) NA else paste0(x,", "'_",i, "')})")
            eval(parse(text=text2))            

            ## Values of the moderators with definition variables
            text3 <- paste0("Sx",i, " <- as.mxMatrix(Sx[[",
                            i, "]], name='Sx",i, "')")
            eval(parse(text=text3))
        }  ## from 'if else' 

        ## Not the final Smatrix as it does not include the error variances
        ## S_no_err= S0+S1*Sx1+S2*Sx2+...
        text4 <- paste0("S", seq_len(length(Sx)), "*Sx",  seq_len(length(Sx)),
                        collapse=" + ")   
        text5 <- paste0("SnoErr <- mxAlgebra(S0 + ", text4, ", name='SnoErr')")
        eval(parse(text=text5))
        
    }   ## from 'else' 


    ## No error variance involved for ALL variables
    invS0 <- mxAlgebra( solve(Iden-Amatrix)%&%SnoErr, name="invS0")

    ## Setup for error variances for mediators and dvs
    ## Excluded level1 as they are IVs
    for (i in 1:length(Level)) {
        ## Select the diagonals of mediators
        text1 <- paste0("F",i, " <- as.mxMatrix(diag(Level[[",i,"]]), name='F",i,"')")
        eval(parse(text=text1))

        ## Extract the error variances of mediators based on the previous
        ## model implied covariance matrix (i-1)
        text2 <- paste0("E",i, " <- mxAlgebra(vec2diag(diag2vec(Iden-invS",i-1,"))*F",
                        i, ", name='E", i, "')")
        eval(parse(text=text2))

        ## Implied covariance matrix with estimated error variances on mediators for (i)
        text3 <- paste("E", 1:i, sep="", collapse="+")
        text4 <- paste("invS",i, " <- mxAlgebra(solve(Iden-Amatrix)%&%(S0+", text3,
                       "), name='invS",i,"')", sep="")
        eval(parse(text=text4))
    }

    ## Final Smatrix including all error variances
    text5 <- paste("E", seq_len(length(Level)), sep="", collapse="+")

    ## Final Smatrix=SnoErr +E1+E2+E3...
    text6 <- paste0("Smatrix <- mxAlgebra(SnoErr+", text5, ", name='Smatrix')")
    eval(parse(text=text6))

    ## Final impliedR
    impliedR <- eval(parse(text="mxAlgebra((Fmatrix%*%solve(Iden-Amatrix))%&%Smatrix, name='impliedR')"))
    ## Vectorized form of impliedR
    vechsR <- mxAlgebra(t(vechs(impliedR)), name="vechsR")

    ## Output the matrices
    text7 <- paste0("invS", seq_len(length(Level)), collapse=",")
    text8 <- paste0("E", seq_len(length(Level)), collapse=",")
    text9 <- paste0("F", seq_len(length(Level)), collapse=",")
    text10 <- paste0("list(Amatrix, Smatrix, Fmatrix, impliedR, vechsR, A0, S0, SnoErr, invS0, Iden,",
                     text7, ",", text8, ",", text9)

    ## Add moderators if there are any
    if (!is.null(Ax)) {
        text10 <- paste0(text10, ",",
                         paste0("A", seq_len(length(Ax)), collapse=","), ",",
                         paste0("Ax", seq_len(length(Ax)), collapse=","))
    }

    if (!is.null(Sx)) {
        text10 <- paste0(text10, ",",
                         paste0("S", seq_len(length(Sx)), collapse=","), ",",
                         paste0("Sx", seq_len(length(Sx)), collapse=","))
    }

    ## Remember to close it    
    text10 <- paste0(text10, ")")
    out <- eval(parse(text=text10))

    ## Add names of the elements in the list
    names(out) <- sapply(out, function(x) x$name)
    out
}

create.Tau2 <- function(RAM, no.var, RE.type=c("Diag", "Symm"),                           
                        Transform=c("expLog", "sqSD"),
                        RE.startvalues=0.05) {
    if (!missing(RAM)) {
        p <- sum(RAM$F)
        no.var <- p*(p-1)/2
    }
    
    ## It has not been tested yet. But it should work for no.var=1.
    if (no.var<1) stop("'no.var' must be at least 1.\n")
    
    RE.type <- match.arg(RE.type)
    Transform <- match.arg(Transform)

    ## Repeat the starting values if there is one value.
    if (length(RE.startvalues)==1) {
        RE.startvalues <- rep(RE.startvalues, no.var)
    }    

    ## Convert the starting values into positive.
    if (any(RE.startvalues <= 0)) {
        RE.startvalues <- abs(RE.startvalues) + 1e-10
        warning("Not all 'RE.startvalues' are positive.\n")
    }

    ## Convert the starting values into log or sqrt scale.
    switch(Transform,
           expLog = RE.startvalues <- log(RE.startvalues),
           sqSD = RE.startvalues <- sqrt(RE.startvalues))

    ## vecTau1: px1 vector of the log (or SD) Tau1
    ## Cor: pxp corr matrix
    ## Tau2: pxp final variance component of variances
    switch(RE.type,
           Symm = {              
               vecTau1 <- create.mxMatrix(paste0(RE.startvalues, "*Tau1_", seq(no.var)),
                                          ncol=1, nrow=no.var, name="vecTau1")
               Cor <- create.mxMatrix(vechs(outer(seq(no.var), seq(no.var),
                                            function(x,y) paste0("0*Cor_", x, "_", y))),
                                      type="Stand", ncol=no.var, nrow=no.var,
                                      lbound=-1, ubound=1, name="Cor")},
           Diag = {
               vecTau1 <- create.mxMatrix(paste0(RE.startvalues, "*Tau1_", seq(no.var)),
                                          ncol=1, nrow=no.var, name="vecTau1")
               Cor <- as.mxMatrix(diag(no.var), name="Cor")})

    switch(Transform,
           expLog = Tau2 <- mxAlgebra( vec2diag(exp(vecTau1)) %&% Cor, name="Tau2" ),
           sqSD = Tau2 <- mxAlgebra( vec2diag(vecTau1 %^% 2) %&% Cor, name="Tau2" ))

    ## Vectorized version of Tau2
    vechTau2 <- mxAlgebra(vech(Tau2), name="vechTau2")
    
    out <- list(Tau2, vecTau1, Cor, vechTau2)

    ## Add names of the elements in the list
    names(out) <- sapply(out, function(x) x$name)
    out
}    

## x: labels of variables by column major
create.V <- function(x, type=c("Symm", "Diag", "Full"), as.mxMatrix=TRUE) {
    type <- match.arg(type)

    switch(type,
           Symm = out <- vec2symMat(x, diag=TRUE, byrow=FALSE),
           Diag = out <- Diag(x),
           Full = {
               no.var <- sqrt(length(x))
               if (no.var%%1 != 0) {
                   stop("Length of '", x, "' is not a square.\n")
               } else {
                   out <- matrix(x, nrow=no.var, ncol=no.var, byrow=FALSE)
               }
               })
           
    out <- apply(out, c(1, 2),
                 function(x) { if (x=="0") x else paste0("0*data.", x) })

    if (as.mxMatrix) out <- as.mxMatrix(out, name="V")
    out
}

## osmasem <- function(model.name="osmasem", Mmatrix, Tmatrix, data, ...) {
##     ## Create known sampling variance covariance matrix
##     Vmatrix <- create.V(data$vlabels, type="Symm", as.mxMatrix=TRUE)
    
##     mxModel(model=model.name, mxData(observed=data$data, type="raw"),
##             mxExpectationNormal(covariance="expCov", means="vechsR",
##                                 dimnames=data$ylabels),
##             mxFitFunctionML(), Mmatrix, Tmatrix, Vmatrix, 
##             mxAlgebra(Tau2+V, name="expCov"),
##             mxCI(c("Amatrix", "Smatrix", "Tau2")), ...)    
## }
    
osmasem <- function(model.name="osmasem", Mmatrix, Tmatrix, data,
                  intervals.type = c("z", "LB"), mxModel.Args=NULL,
                  mxRun.Args=NULL, suppressWarnings=TRUE,
                  silent=FALSE, run=TRUE, ...) {

    intervals.type <- match.arg(intervals.type)
    switch(intervals.type,
           z = intervals <- FALSE,
           LB = intervals <- TRUE)
    
    ## Create known sampling variance covariance matrix
    Vmatrix <- create.V(data$vlabels, type="Symm", as.mxMatrix=TRUE)
    
    mx.model <- eval(parse(text="mxModel(model=model.name,
                        mxData(observed=data$data, type='raw'),
                        mxExpectationNormal(covariance='expCov',
                                            means='vechsR',
                                            dimnames=data$ylabels),
                        mxFitFunctionML(), Mmatrix, Tmatrix, Vmatrix,
                        mxAlgebra(Tau2+V, name='expCov'),
                        mxCI(c('Amatrix', 'Smatrix', 'Tau2')))"))

    ## Add additional arguments to mxModel
    if (!is.null(mxModel.Args)) {
        for (i in seq_along(mxModel.Args)) {
            mx.model <- mxModel(mx.model, mxModel.Args[[i]])
        }
    }

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

    out <- list(call=match.call(), Mmatrix=Mmatrix, Tmatrix=Tmatrix,
                Vmatrix=Vmatrix, data=data, mxModel.Args=mxModel.Args,
                mxRun.Args=mxRun.Args, mx.model=mx.model, mx.fit=mx.fit)
    class(out) <- 'osmasem'
    out
}

coef.osmasem <- function(object, select=c("fixed", "all", "random"), ...) {
    if (!is.element("osmasem", class(object)))
        stop("\"object\" must be an object of class \"osmasem\".")
    
    mx.coef <- coef(object$mx.fit)
    ## index for untransformed random effects (not the correct ones!) 
    index <- grep("Tau1_|Cor_", names(mx.coef))

    select <- match.arg(select)
    switch(select,
           fixed =  mx.coef <- mx.coef[-index],
           random = mx.coef <- mx.coef[index])
    if (select!="fixed") warning("\"Tau1_xx\" is not the variance component of the random effects. Please use VarCorr().")    
    mx.coef    
}

vcov.osmasem <- function(object, select=c("fixed", "all", "random"), ...) {
    if (!is.element("osmasem", class(object)))
        stop("\"object\" must be an object of class \"osmasem\".")
    
    mx.vcov <- vcov(object$mx.fit)
    ## index for untransformed random effects (not the correct ones!) 
    index <- grep("Tau1_|Cor_", rownames(mx.vcov))

    select <- match.arg(select)
    switch(select,
           fixed =  mx.vcov <- mx.vcov[-index, -index],
           random = mx.vcov <- mx.vcov[index, index])
    if (select!="fixed") warning("\"Tau1_xx\" is not the variance component of the random effects.")
    mx.vcov 
}

VarCorr <- function(x, ...) {
    if (!(class(x) %in% c("meta", "osmasem")))
        stop("\"x\" must be an object of either class \"meta\" or \"osmasem\".")

    if (class(x)=="meta") {
        eval(parse(text="mxEval(Tau, x$mx.fit)"))
    } else {
        out <- eval(parse(text="mxEval(Tau2, x$mx.fit)"))

        ## No. of variables in the model
        p <- ncol(out)
        names.tau2 <- paste0("Tau2_", seq_len(p))    
        dimnames(out) <- list(names.tau2, names.tau2)
        out
    }     
}


## Fit a saturated model with either a diagonal or symmetric variance component of random effects
.osmasemSatIndMod <- function(osmasem.obj=NULL, model=c("Saturated", "Independence"),
                              Std.Error=FALSE, Tmatrix, data, mxModel.Args=NULL,
                              mxRun.Args=NULL, suppressWarnings=TRUE,
                              silent=FALSE, run=TRUE, ...) {

    if (!is.null(osmasem.obj)) {
        if (!is.element("osmasem", class(osmasem.obj)))
            stop("\"osmasem.obj\" must be an object of class \"osmasem\".")

        Tmatrix <- osmasem.obj$Tmatrix
        Vmatrix <- osmasem.obj$Vmatrix
        ## No. of variables in the model
        p <- nrow(Tmatrix$Cor$values)
        data <- osmasem.obj$data
    } else {
        p <- length(data$ylabels)
        ## Create known sampling variance covariance matrix
        Vmatrix <- create.V(data$vlabels, type="Symm", as.mxMatrix=TRUE)
    }    

    model <- match.arg(model)

    if (model=="Saturated") {
        ## Saturated model
        Mmatrix <- mxMatrix(type="Full", free=TRUE, labels=paste0("Mu", seq_len(p)),
                            nrow=1, ncol=p, name="Mu")
    } else {
        ## Independence model
        Mmatrix <- mxMatrix(type="Full", free=FALSE, labels=paste0("Mu", seq_len(p)),
                            nrow=1, ncol=p, name="Mu")        
    }
    
    mx.model <- eval(parse(text="mxModel(model=model,
                        mxData(observed=data$data, type='raw'),
                        mxExpectationNormal(covariance='expCov',
                                            means='Mu',
                                            dimnames=data$ylabels),
                        mxFitFunctionML(), Mmatrix, Tmatrix, Vmatrix,
                        mxAlgebra(Tau2+V, name='expCov'))"))

    ## Add additional arguments to mxModel
    if (!is.null(mxModel.Args)) {
        for (i in seq_along(mxModel.Args)) {
            mx.model <- mxModel(mx.model, mxModel.Args[[i]])
        }
    }

    if (Std.Error==FALSE) {
        mx.model <- mxOption(mx.model, "Calculate Hessian", "No") 
        mx.model <- mxOption(mx.model, "Standard Errors"  , "No")
    }
    
    ## Return mx model without running the analysis
    if (run==FALSE) {
        return(mx.model)
    } else {        
        mx.fit <- tryCatch(do.call(mxRun, c(list(mx.model, silent=silent,
                                                 suppressWarnings=suppressWarnings),
                                            mxRun.Args)), error = function(e) e)
    }
    
    # try to run it with error message as output
    if (inherits(mx.fit, "error")) {
        warning(print(mx.fit))
    }
    mx.fit
}

anova.osmasem <- function(object, ..., all=FALSE) {
  base <- lapply(list(object), function(x) x$mx.fit)
  comparison <- lapply(list(...), function(x) x$mx.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

summary.osmasem <- function(object, Saturated=FALSE, numObs, ...) {
    if (!is.element("osmasem", class(object)))
        stop("\"object\" must be an object of class \"osmasem\".")

    ## Ind.stat <- summary(.osmasemSatIndMod(object, model="Independence",
    ##                                    Std.Error=FALSE, silent=TRUE))   

    ## If numObs is not provided, use total N
    if (missing(numObs)) {
        ## numObs <- sum(object$data$n)-length(object$data$n)
        numObs <- sum(object$data$n)
    }

    ## Calculate chi-square statistic
    if (Saturated) {
        Sat.stat <- summary(.osmasemSatIndMod(object, model="Saturated",
                                              Std.Error=FALSE, silent=TRUE))
        out <- summary(object$mx.fit,
                       SaturatedLikelihood=Sat.stat$Minus2LogLikelihood,
                       SaturatedDoF=Sat.stat$degreesOfFreedom,
                       numObs=numObs, ...)        
    } else {
        out <- summary(object$mx.fit,
                       numObs=numObs, ...)
    }
    
    out$parameters$`z value` <- with(out$parameters, Estimate/Std.Error)
    out$parameters$`Pr(>|z|)` <- 2*(1-pnorm(abs(out$parameters$`z value`)))
    out
}

osmasemR2 <- function(object, R2.truncate=TRUE) {
    if (!is.element("osmasem", class(object)))
        stop("\"object\" must be an object of class \"osmasem\".")

    ## No. of variables in the model
    p <- nrow(object$Tmatrix$Cor$values)

    outer(seq_len(p), seq_len(p), function(x, y) paste0("Tau2_", x, "_", y))
    
    Sat.model <- .osmasemSatIndMod(object, Std.Error=FALSE)
    Tau2.0 <- diag(eval(parse(text = "mxEval(Tau2, Sat.model)")))       
    Tau2.1 <- diag(eval(parse(text = "mxEval(Tau2, object$mx.fit)")))

    R2 <- (Tau2.0-Tau2.1)/Tau2.0

    ## Truncate negative R2 to 0
    if (R2.truncate) R2[R2<0] <- 0
        
    names.tau2 <- paste0("Tau2_", seq_len(p), "_", seq_len(p))
    names(R2) <- names(Tau2.0) <- names(Tau2.1) <- names.tau2
    
    list(Tau2.0=Tau2.0, Tau2.1=Tau2.1, R2=R2)
}

## RAMmodelV <- function(Amatrix, Smatrix, Fmatrix, Mmatrix, Vmatrix,
##                       data, obs.var, ...) {

##     Sfull <- mxAlgebra(Smatrix + Vmatrix, name="Sfull")
##     fitfun <- mxExpectationRAM(A="Amatrix", S="Sfull", F="Fmatrix",
##                                M="Mmatrix", dimnames=obs.var)
##     fit <- mxFitFunctionML()
##     mxModel(mxData(observed=data, type="raw"),
##             Amatrix, Smatrix, Fmatrix, Mmatrix, Vmatrix, Sfull,
##             fitfun, fit, ...)
## }

## NORMmodelV <- function(Mu, Sigma, Fmatrix, Vmatrix, data, obs.var, ...) {

##     Sfull <- mxAlgebra(Sigma + Vmatrix, name="Sfull")
##     fitfun <- mxExpectationNormal(covariance="Sfull",
##                                   means=Mu@name,
##                                   dimnames=obs.var)
##     fit <- mxFitFunctionML()
##     mxModel(mxData(observed=data, type="raw"),
##             Mu, Sigma, Fmatrix, Vmatrix, Sfull, fitfun, fit, ...)
## }

