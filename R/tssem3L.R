tssem3L1 <- function(Cov, n, cluster, RE.typeB=c("Diag", "Symm"),
                     RE.typeW=c("Diag", "Symm"),
                     acov=c("weighted", "individual", "unweighted"), ...) {

    if (length(Cov)!=length(n)) stop("\"Cov\" and \"n\" are in different lengths!\n")
    if (length(Cov)!=length(cluster)) stop("\"Cov\" and \"cluster\" are in different lengths!\n")
    
    RE.typeB <- match.arg(RE.typeB, c("Diag", "Symm"))
    RE.typeW <- match.arg(RE.typeW, c("Diag", "Symm"))    
    acov <- match.arg(acov, c("weighted", "individual", "unweighted"))

    ## Save the original cluster (data)
    cluster_data <- cluster
    ## Cluster variable name
    cluster <- "cluster"
    
    df <- Cor2DataFrame(x=Cov, n=n, acov=acov)
    ## Add cluster to the dataset
    df$data <- cbind(df$data, cluster=as.factor(cluster_data))

    ## No. of observed variables in RAM
    p <- nrow(Cov[[1]])
    
    ## Dimension of random effects: no.var by no.var
    ylabels <- df$ylabels
    p_star <- length(ylabels)
    
    ## Saturated model
    impliedR <- mxMatrix(type="Stand", nrow=p, ncol=p, free=TRUE,
                         values=unname(apply(df$data[, ylabels], 2, mean, na.rm=TRUE)),
                         labels=ylabels, name="impliedR")

    ## dimnames are required
    vechsR <- mxAlgebra(t(vechs(impliedR)), dimnames=list(NULL, ylabels), name="vechsR")    
    Mmatrix <- list(vechsR=vechsR, impliedR=impliedR)
        
    TmatrixB <- create.Tau2(no.var=p_star, RE.type=RE.typeB, Transform="expLog", level="between")
    TmatrixW <- create.Tau2(no.var=p_star, RE.type=RE.typeW, Transform="expLog", level="within")

    mx.fit <- osmasem3L(cluster="cluster", data=df, Mmatrix=Mmatrix, TmatrixB=TmatrixB,
                        TmatrixW=TmatrixW)
    class(mx.fit) <- c(class(mx.fit), "tssem3L1")
    mx.fit
}

tssem3L2 <- function(tssem1.obj, RAM=NULL, Amatrix=NULL, Smatrix=NULL, Fmatrix=NULL,
                     diag.constraints=FALSE, intervals.type = c("z", "LB"),
                     mx.algebras=NULL, mxModel.Args=NULL, model.name=NULL,
                     suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) {
    if ( !is.element( class(tssem1.obj)[2], "tssem3L1") )
        stop("\"tssem1.obj\" must be of class \"tssem3L1\".")

    Cov <- vec2symMat(coef(tssem1.obj, select="fixed"), diag=FALSE)
    dimnames(Cov) <- list(tssem1.obj$data$obslabels, tssem1.obj$data$obslabels)

    wls(Cov=Cov, aCov=vcov(tssem1.obj), n=sum(tssem1.obj$data$n),
        RAM=RAM, Amatrix=Amatrix, Smatrix=Smatrix, Fmatrix=Fmatrix,
        diag.constraints=diag.constraints, cor.analysis=TRUE,
        intervals.type=intervals.type, mx.algebras=mx.algebras, mxModel.Args=mxModel.Args,
        model.name=model.name, suppressWarnings=suppressWarnings,
        silent=silent, run=run, ...) 
}

    

