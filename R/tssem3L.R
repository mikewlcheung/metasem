tssem3L1 <- function(cluster=NULL, RE.typeB=c("Diag", "Symm"), RE.typeW=c("Diag", "Symm"),
                     data, ...) {
  
    RE.typeB <- match.arg(RE.typeB, c("Diag", "Symm"))
    RE.typeW <- match.arg(RE.typeW, c("Diag", "Symm"))    
 
    ## No. of observed variables in RAM
    p <- length(data$obslabels)
    
    ## Dimension of random effects: no.var by no.var
    ylabels <- data$ylabels
    p_star <- length(ylabels)
    
    ## Saturated model
    impliedR <- mxMatrix(type="Stand", nrow=p, ncol=p, free=TRUE,
                         values=unname(apply(data$data[, ylabels, drop=FALSE], 2, mean, na.rm=TRUE)),
                         labels=ylabels, name="impliedR")

    ## dimnames are required
    vechsR <- mxAlgebra(t(vechs(impliedR)), dimnames=list(NULL, ylabels), name="vechsR")    
    Mmatrix <- list(vechsR=vechsR, impliedR=impliedR)
        
    TmatrixB <- create.Tau2(no.var=p_star, RE.type=RE.typeB, Transform="expLog", level="between")
    TmatrixW <- create.Tau2(no.var=p_star, RE.type=RE.typeW, Transform="expLog", level="within")

    ## RE.typeW is required in osmasem3L
    mx.fit <- osmasem3L(model.name="tssem3L1", cluster=cluster, data=data, Mmatrix=Mmatrix, TmatrixB=TmatrixB,
                        TmatrixW=TmatrixW, RE.typeW=RE.typeW)
    class(mx.fit) <- c(class(mx.fit), "tssem3L1")
    mx.fit
}

tssem3L2 <- function(tssem1.obj, RAM=NULL, Amatrix=NULL, Smatrix=NULL, Fmatrix=NULL,
                     diag.constraints=FALSE, intervals.type = c("z", "LB"),
                     mx.algebras=NULL, mxModel.Args=NULL, model.name=NULL,
                     suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {
    if ( !is.element( class(tssem1.obj)[2], "tssem3L1") )
        stop("\"tssem1.obj\" must be of class \"tssem3L1\".")

    # Cov <- vec2symMat(coef(tssem1.obj, select="fixed"), diag=FALSE)
    Cov <- tssem1.obj$mx.fit$impliedR$values
    dimnames(Cov) <- list(tssem1.obj$data$obslabels, tssem1.obj$data$obslabels)

    wls(Cov=Cov, aCov=vcov(tssem1.obj), n=sum(tssem1.obj$data$n),
        RAM=RAM, Amatrix=Amatrix, Smatrix=Smatrix, Fmatrix=Fmatrix,
        diag.constraints=diag.constraints, cor.analysis=TRUE,
        intervals.type=intervals.type, mx.algebras=mx.algebras, mxModel.Args=mxModel.Args,
        model.name=model.name, suppressWarnings=suppressWarnings,
        silent=silent, run=run, ...) 
}

    

