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
create.vechsR <- function(A0, S0, F0=NULL, Ax=NULL, Sx=NULL,
                          A.lbound=NULL, A.ubound=NULL) {

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

            ## Fix the diagonals of S1, S2... to 0
            text4 <- paste0("diag(S",i,"$values) <- 0")
            eval(parse(text=text4))
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

    ## Add lbound and ubound to A0
    if (!is.null(A.lbound)) {
        if (is.matrix(A.lbound)) {
            if (dim(A.lbound)==dim(out$A0$values)) {
                out$A0$lbound <- A.lbound
            } else {
                stop("The dimensions of 'A.lbound' do not match those in 'A0'.\n")
            }
        } else {  # lbound is not matrix
            out$A0$lbound  <- matrix(A.lbound, nrow=nrow(out$A0$values),
                                     ncol=ncol(out$A0$values),
                                     dimnames=dimnames(out$A0$values))
        }
    }

    if (!is.null(A.ubound)) {
        if (is.matrix(A.ubound)) {
            if (dim(A.ubound)==dim(out$A0$values)) {
                out$A0$ubound <- A.ubound
            } else {
                stop("The dimensions of 'A.ubound' do not match those in 'A0'.\n")
            }
        } else {  # ubound is not matrix
            out$A0$ubound  <- matrix(A.ubound, nrow=nrow(out$A0$values),
                                     ncol=ncol(out$A0$values),
                                     dimnames=dimnames(out$A0$values))
        } 
    }      
                
    out
}

create.Tau2 <- function(RAM, no.var, Tau1.labels=seq(no.var),
                        RE.type=c("Diag", "Symm", "Zero", "User"),
                        level=c("single", "between", "within"),
                        RE.User=NULL, Transform=c("expLog", "sqSD"),
                        RE.startvalues=0.05) {

    level <- match.arg(level)
    switch(level,
           single = level_suffix <- "",
           between= level_suffix <- "B",
           within = level_suffix <- "W")

    if (!missing(RAM)) {
        p <- sum(RAM$F)
        ## no. of correlation coefficients
        no.var <- p*(p-1)/2
    }
    
    ## It has not been tested yet. But it should work for no.var=1.
    if (no.var<1) stop("'no.var' must be at least 1.\n")

    if (length(Tau1.labels) != no.var) stop("Length of \"Tau1.labels\" is different from \"no.var\".\n")
    Tau1.labels <- paste0(Tau1.labels, level_suffix)
    
    RE.type <- match.arg(RE.type)
    Transform <- match.arg(Transform)

    ## Repeat the starting values if there is only one value.
    if (length(RE.startvalues)==1) {
        RE.startvalues <- rep(RE.startvalues, no.var)
    }    

    ## Convert the starting values into positive.
    if (any(RE.startvalues < 0)) {
        RE.startvalues <- abs(RE.startvalues) + 1e-10
        warning("Not all 'RE.startvalues' are positive.\n")
    }

    ## Set the starting values at 0 for RE.type="Zero"
    if (RE.type=="Zero") RE.startvalues <- rep(0, no.var)
    
    ## Convert the starting values into log or sqrt scale.
    switch(Transform,
           expLog = RE.startvalues <- log(RE.startvalues),
           sqSD = RE.startvalues <- sqrt(RE.startvalues))

    ## vecTau1: px1 vector of the log (or SD) Tau1
    ## Cor: pxp corr matrix
    ## Tau2: pxp final variance component of variances
    switch(RE.type,
           Symm = {              
               vecTau1 <- create.mxMatrix(paste0(RE.startvalues, "*Tau1_", Tau1.labels),
                                          ncol=1, nrow=no.var, name="vecTau1")
               Cor <- create.mxMatrix(vechs(outer(seq(no.var), seq(no.var),
                                            function(x,y) paste0("0*Cor_", x, "_", y, level_suffix))),
                                      type="Stand", ncol=no.var, nrow=no.var,
                                      lbound=-0.99, ubound=0.99, name=paste0("Cor", level_suffix))},
           Diag = {
               vecTau1 <- create.mxMatrix(paste0(RE.startvalues, "*Tau1_", Tau1.labels),
                                          ncol=1, nrow=no.var, name="vecTau1")
               ## Uncorrelated
               Cor <- as.mxMatrix(diag(no.var), name=paste0("Cor", level_suffix))},
           Zero = {
               vecTau1 <- create.mxMatrix(RE.startvalues, type="Full", ncol=1,
                                          nrow=no.var, name="vecTau1")
               Cor <- as.mxMatrix(diag(no.var), name=paste0("Cor", level_suffix))},
           User = { if (is.null(RE.User)) stop("'RE.User', the ", no.var, " by ", no.var,
                                               " symmetric matrix of TRUE or FALSE on the variance componment, must be specified.\n")
                                               ## Check the symmetry of RE.User
               if (!isSymmetric(RE.User)) stop("'RE.User' is not symmetric.\n")
               if (!all(dim(RE.User)==c(no.var, no.var))) stop("The dimensions of 'RE.User' are different 'no.var'.\n")
               ## Check if the covariances are free but the variances are fixed.                                
               if (no.var>1) {
                   for (i in 1:(no.var-1)) {
                       for (j in (i+1):no.var) {
                           if (RE.User[i,j] & !all(RE.User[i,i], RE.User[j,j])) {
                               stop("RE.User[",j,",",i,"] is free but either RE.User[",i,",",
                                    i,"] or RE.User[",j,",",j,"] is fixed.\n") }}}}
                                               
               vecTau1 <- paste0(RE.startvalues, "*Tau1_", Tau1.labels)
               ## Fixed a bug that RE.startvalues are not used when the variances are not free.
               ## vecTau1[diag(RE.User)==FALSE] <- 0                                
               switch(Transform,
                      expLog = vecTau1[diag(RE.User)==FALSE] <- log(0),
                      sqSD =   vecTau1[diag(RE.User)==FALSE] <- sqrt(0))
               vecTau1 <- create.mxMatrix(vecTau1, ncol=1, nrow=no.var, name="vecTau1")
               Cor <- outer(seq(no.var), seq(no.var),
                            function(x,y) paste0("0*Cor_", x, "_", y, level_suffix))
               Cor[RE.User==FALSE] <- 0
               Cor <- create.mxMatrix(vechs(Cor), type="Stand", ncol=no.var, nrow=no.var,
                                      lbound=-0.99, ubound=0.99, name=paste0("Cor", level_suffix))}
               )

    switch(Transform,
           expLog = Tau2 <- eval(parse(text=paste0("mxAlgebra( vec2diag(sqrt(exp(vecTau1))) %&% ",
                                                   "Cor", level_suffix,
                                                   ", name='Tau2", level_suffix, "')"))) ,
           sqSD = Tau2 <- eval(parse(text=paste0("mxAlgebra( vec2diag(vecTau1) %&% ",
                                                   "Cor", level_suffix,
                                                   ", name='Tau2", level_suffix, "')"))))

    ## Vectorized version of Tau2
    vechTau2 <- eval(parse(text=paste0("mxAlgebra(vech(Tau2", level_suffix,
                                       "), name='vechTau2", level_suffix, "')")))
    
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
    
osmasem <- function(model.name="osmasem", RAM=NULL, Mmatrix=NULL,
                    Tmatrix=NULL, Jmatrix=NULL, Ax=NULL, Sx=NULL,
                    A.lbound=NULL, A.ubound=NULL,
                    RE.type=c("Diag", "Symm"), data,
                    subset.variables=NULL, subset.rows=NULL, 
                    intervals.type = c("z", "LB"),
                    mxModel.Args=NULL, mxRun.Args=NULL,
                    suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {

    RE.type <- match.arg(RE.type)
    
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

    ## If RAM is provided, create the matrices based on it.
    if (!is.null(RAM)) {
        Mmatrix <- create.vechsR(A0=RAM$A, S0=RAM$S, F0=RAM$F, Ax=Ax, Sx=Sx,
                                 A.lbound=A.lbound, A.ubound=A.ubound)
        Tmatrix <- create.Tau2(RAM=RAM, RE.type=RE.type, Transform="expLog")
    }

    ## Create known sampling variance covariance matrix
    Vmatrix <- create.V(vlabels, type="Symm", as.mxMatrix=TRUE)

    ## If Jmatrix is not specified, create an identify matrix
    ## Use Vmatrix$Cor$values as it is an identity matrix
    if (is.null(Jmatrix)) {
        Jmatrix <- Vmatrix$values
        diag(Jmatrix) <- 1
        Jmatrix <- as.mxMatrix(Jmatrix, name="Jmatrix")
    } else {
        Jmatrix$name  <- "Jmatrix"
    }
    
    if (is.null(subset.rows)) {
        subset.rows <- rep(TRUE, nrow(data$data))
    }
    ## Select data for the analaysis
    mx.data <- data$data[subset.rows, ]

    ## Dirty trick to avoid the warning of "no visible binding for global variable"
    mx.model <- eval(parse(text="mxModel(model=model.name,
                           mxData(observed=mx.data, type='raw'),
                           mxExpectationNormal(covariance='expCov',
                                               means='vechsR',
                                               dimnames=ylabels),
                           mxFitFunctionML(),
                           Mmatrix, Tmatrix, Vmatrix, Jmatrix,
                           mxAlgebra(Jmatrix %&% Tau2 + V, name='expCov'),
                           mxCI(c('Amatrix', 'Smatrix', 'Tau2')))"))
    
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

    ## Assign same starting values, e.g., A and S matrices are different from those in provided in Jmatrix
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
                Tmatrix=Tmatrix,
                Vmatrix=Vmatrix, 
                Jmatrix=Jmatrix, 
                data=data,
                labels=list(obslabels=obslabels, ylabels=ylabels, vlabels=vlabels),
                mxModel.Args=mxModel.Args,
                mxRun.Args=mxRun.Args,
                subset.variables=subset.variables,
                subset.rows=subset.rows,
                mx.model=mx.model, mx.fit=mx.fit)
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
    if (select!="fixed") warning("\"Tau1_xx\" is not the variance component of the random effects.\nPlease use VarCorr() to get the variance component.\n")    
    mx.coef    
}

vcov.osmasem <- function(object, select=c("fixed", "all", "random"),
                         robust=FALSE, ...) {
    if (!is.element("osmasem", class(object)))
        stop("\"object\" must be an object of class \"osmasem\".")

    # labels of the parameters    
    ## my.name <- summary(object$mx.fit)$parameters$name
    my.name <- names(omxGetParameters(object$mx.fit))
    my.name <- my.name[!is.na(my.name)]

    ## index for untransformed random effects (not the correct ones!) 
    index.random <- grep("Tau1_|Cor_", my.name) 
    
    select <- match.arg(select)
    switch( select,
         ## all = my.name <- my.name,
         fixed =  my.name <- my.name[-index.random],
         random = my.name <- my.name[index.random]
         )

    if (robust) {
        out <- suppressMessages(imxRobustSE(object$mx.fit, details=TRUE))$cov
    } else {
        out <- vcov(object$mx.fit)
    }
    if (select!="fixed") warning("\"Tau1_xx\" is not the variance component of the random effects.")
    out[my.name, my.name]
}

VarCorr <- function(x, ...) {
    if (all(!is.element(c("meta", "osmasem", "osmasem3L"), class(x))))
        stop("\"x\" must be an object of either class \"meta\", \"osmasem\" or \"osmasem3L\".")

    switch(class(x)[1],
           meta = out <- eval(parse(text="mxEval(Tau, x$mx.fit)")),
           osmasem = {
               out <- eval(parse(text="mxEval(Tau2, x$mx.fit)"))

               ## ## No. of variables in the model
               ## p <- ncol(out)
               ## names.tau2 <- paste0("Tau2_", seq_len(p))

               ## Get the names Tau1_1, Tau1_2... or Tau1_a, Tau1_b... from x$Tmatrix
               ## Replace Tau1 with Tau2
               names.tau2 <- gsub("Tau1", "Tau2", c(x$Tmatrix$vecTau1$labels))
               dimnames(out) <- list(names.tau2, names.tau2)},
           osmasem3L = {
               ## Tau2 (within)
               Tau2W <- eval(parse(text="mxEval(Tau2W, x$mx.fit)"))
               ## Labels ended with "W"
               names.tau2W <- gsub("W$", "", c(x$TmatrixW$vecTau1$labels))
               names.tau2W <- gsub("Tau1", "Tau2", names.tau2W)               
               dimnames(Tau2W) <- list(names.tau2W, names.tau2W)

               ## Tau2 (between)
               Tau2B <- eval(parse(text="mxEval(B.Tau2B, x$mx.fit)"))
               ## Labels ended with "B" 
               names.tau2B <- gsub("B$", "", c(x$TmatrixB$vecTau1$labels))
               names.tau2B <- gsub("Tau1", "Tau2", names.tau2B)
               dimnames(Tau2B) <- list(names.tau2B, names.tau2B)
               out <- list(Tau2W=Tau2W, Tau2B=Tau2B)},
           stop("Invalid class:", class(x)))

    out    
}


## Fit a saturated model with either a diagonal or symmetric variance component of random effects
.osmasemSatIndMod <- function(osmasem.obj=NULL, model=c("Saturated", "Independence"),
                              Std.Error=FALSE, Tmatrix, data, 
                              mxModel.Args=NULL, mxRun.Args=NULL, suppressWarnings=TRUE,
                              silent=FALSE, run=TRUE, ...) {

    if (!is.null(osmasem.obj)) {
        if (!is.element("osmasem", class(osmasem.obj)))
            stop("\"osmasem.obj\" must be an object of class \"osmasem\".")

        if (missing(Tmatrix)) {
            Tmatrix <- osmasem.obj$Tmatrix
        }                
        ## No. of variables in the model
        p <- nrow(Tmatrix$Cor$values)
        data <- osmasem.obj$data
        ylabels <- osmasem.obj$labels$ylabels
        Vmatrix <- osmasem.obj$Vmatrix        
    } else {
        ## FIXME: It does not work.
        ## If there is no osmasem.obj, create a new saturated/independence model from data.
        ylabels <- data$ylabels
        p <- length(ylabels)
        ## Create known sampling variance covariance matrix
        Vmatrix <- create.V(data$labels$vlabels, type="Symm", as.mxMatrix=TRUE)
    }    

    model <- match.arg(model)

    if (model=="Saturated") {
        ## Saturated model
        st.values <- apply(data$data[, ylabels, drop=FALSE], 2,  mean, na.rm=TRUE)
        Mmatrix <- mxMatrix(type="Full", free=TRUE, values=st.values,
                            labels=paste0("Mu", seq_len(p)),
                            nrow=1, ncol=p, dimnames=list(NULL, ylabels), name="Mu")
        
        Mmatrix <- mxMatrix(type="Full", free=TRUE, labels=paste0("Mu", seq_len(p)),
                            nrow=1, ncol=p, name="Mu")
    } else {
        ## Independence model
        Mmatrix <- mxMatrix(type="Full", free=FALSE, labels=paste0("Mu", seq_len(p)),
                            nrow=1, ncol=p, name="Mu")        
    }

    ## Filter the data for the analysis
    mx.data <- data$data[osmasem.obj$subset.rows, ]

    ## Dirty trick to avoid the warning of "no visible binding for global variable"
    mx.model <- eval(parse(text="mxModel(model=model,
                                 mxData(observed=mx.data, type='raw'),
                                 mxExpectationNormal(covariance='expCov',
                                                     means='Mu',
                                                     dimnames=osmasem.obj$labels$ylabels),
                                 mxFitFunctionML(), Mmatrix, Tmatrix, Vmatrix,
                                 mxAlgebra(Tau2+V, name='expCov'))"))
    
    ## Add additional arguments to mxModel
    if (!is.null(mxModel.Args)) {
        for (i in seq_along(mxModel.Args)) {
            mx.model <- mxModel(mx.model, mxModel.Args[[i]])
        }
    }

    ## It may speed up the analysis
    if (Std.Error==FALSE) {
        mx.model <- mxOption(mx.model, "Calculate Hessian", "No")
        mx.model <- mxOption(mx.model, "Standard Errors", "No")
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
    ## if (inherits(mx.fit, "error")) {
    ##     warning(print(mx.fit))
    ## }
    mx.fit
}

anova.osmasem <- function(object, ..., all=FALSE) {
  base <- lapply(list(object), function(x) x$mx.fit)
  comparison <- lapply(list(...), function(x) x$mx.fit)
  mxCompare(base=base, comparison=comparison, all=all)
}

summary.osmasem <- function(object, fitIndices=FALSE, numObs,
                            robust=FALSE, ...) {
    if (!is.element("osmasem", class(object)))
        stop("\"object\" must be an object of class \"osmasem\".")

    ## Ind.stat <- summary(.osmasemSatIndMod(object, model="Independence",
    ##                                    Std.Error=FALSE, silent=TRUE))   

    ## If numObs is not provided, use the total N
    if (missing(numObs)) {
        ## numObs <- sum(object$data$n)-length(object$data$n)
        numObs <- sum(object$data$n[object$subset.rows])
    }

    ## Calculate chi-square statistic and other fit indices
    if (fitIndices)
    {   ## Check if the dimensions of Tmatrix and Vmatrix match.
        ## They can be different for rcmasem
        ## No. of effect sizes
        p <- dim(object$Vmatrix$values)[1]
        
        if (length(object$Tmatrix$vecTau1$labels) != p)
        {   ## Check if there are correlations,
            ## If yes -> Symm, no -> Diag
            ## if (sum(object$Tmatrix$Cor$free)>0) RE.type="Symm" else RE.type="Diag"
            ## T0 <- create.Tau2(no.var=p, RE.type=RE.type)

            ## TODO: For rcmasem, the variance component cannot be Diag;
            ## otherwise, the saturated model fits poorer than the rcmasem.
            T0 <- create.Tau2(no.var=p, RE.type="Symm")
            
            Sat.stat <- .osmasemSatIndMod(object, model="Saturated", Tmatrix=T0,
                                          Std.Error=FALSE, silent=TRUE)
            Ind.stat <- .osmasemSatIndMod(object, model="Independence", Tmatrix=T0,
                                          Std.Error=FALSE, silent=TRUE)            
        } else {
            Sat.stat <- .osmasemSatIndMod(object, model="Saturated",
                                          Std.Error=FALSE, silent=TRUE)
            Ind.stat <- .osmasemSatIndMod(object, model="Independence", 
                                          Std.Error=FALSE, silent=TRUE) 
        }
       
        ## Sat.model=FALSE for saturated model if there are either errors or nonconvergent.
        if ( inherits(Sat.stat, "error") | !(Sat.stat$output$status$code %in% c(0,1)) )
        {
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
        if ( inherits(Ind.stat, "error") |
             !(Ind.stat$output$status$code %in% c(NA, 0, 1)) ) 
        {
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
        if (Sat.model)
        {
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

    ## Modified the SE with robust SE
    if (robust) {
        robust.SE <- suppressMessages(imxRobustSE(object$mx.fit))
        out$parameters[, "Std.Error"] <- robust.SE[out$parameters$name]
    }

    ## Additional output to the mx summary
    out$parameters$`z value` <- with(out$parameters, Estimate/Std.Error)
    out$parameters$`Pr(>|z|)` <- 2*(1-pnorm(abs(out$parameters$`z value`)))
    out
}

osmasemR2 <- function(model1, model0, R2.truncate=TRUE) {
    if (!all(c(class(model0), class(model1)) %in% "osmasem"))
        stop("Both \"model0\" and \"model1\" must be objects of class \"osmasem\".")

    ## Tau2.0 <- diag(eval(parse(text = "mxEval(Tau2, model0$mx.fit)")))       
    ## Tau2.1 <- diag(eval(parse(text = "mxEval(Tau2, model1$mx.fit)")))
    Tau2.0 <- diag(VarCorr(model0))
    Tau2.1 <- diag(VarCorr(model1))               

    R2 <- (Tau2.0-Tau2.1)/Tau2.0

    ## Truncate negative R2 to 0
    if (R2.truncate) R2[R2<0] <- 0
    
    p <- nrow(model0$Tmatrix$Cor$values)        
    names.tau2 <- paste0("Tau2_", seq_len(p), "_", seq_len(p))
    names(R2) <- names(Tau2.0) <- names(Tau2.1) <- names.tau2
    
    list(Tau2.0=Tau2.0, Tau2.1=Tau2.1, R2=R2)
}

osmasemSRMR <- function(x) {
    ## Check if there are moderators in either A1 or S1
    if (!is.null(x$Mmatrix$A1) | !is.null(x$Mmatrix$S1))
        stop("Moderators are not allowed in calculating the SRMR in OSMASEM.\n")

    ## Saturated model which should be very close to the stage 1 results in the TSSEM
    Sat.stat <- .osmasemSatIndMod(x, model="Saturated", Std.Error=FALSE, silent=TRUE)

    ## Similar to the sample correlation matrix
    sampleR <- vec2symMat(eval(parse(text = "mxEval(Mu, Sat.stat)")), diag=FALSE)

    ## Model implied correlation matrix
    impliedR <- mxEval(impliedR, x$mx.fit)

    ## SRMR (Cheung, 2015 book, Eq. 2.36)
    sqrt(mean(vechs(sampleR-impliedR)^2))
}

## Generate names for OSMASEM
## x: a vector of observed variables
## output: a vector of correlation coefficients and sampling covariance matrix
## This function is used to select data in data created by Cor2DataFrame()
.genCorNames <- function(x) {
  ## Names for the correlation elements
  psNames <- vechs(outer(x, x, paste, sep = "_"))
  
  ## Names for the sampling covariance matrix of the correlation vector
  psCovNames <- paste("C(", outer(psNames, psNames, paste, sep = " "), ")", sep="")
  psCovNames <- vech(matrix(psCovNames, nrow=length(psNames), ncol=length(psNames)))

  list(ylabels=psNames, vlabels=psCovNames)
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
