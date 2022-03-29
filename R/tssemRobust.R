tssemRobust1 <- function(Cov, n, cluster, RE.type=c("Diag", "Symm", "Zero"),
                         acov=c("weighted", "individual", "unweighted"),
                         method=c("ML", "REML"), ...) {
    
    if (!requireNamespace("metafor", quietly=TRUE))    
        stop("\"metafor\" package is required for this function.")

    RE.type <- match.arg(RE.type, c("Diag", "Symm", "Zero"))
    acov <- match.arg(acov, c("weighted", "individual", "unweighted"))
    method <- match.arg(method, c("ML", "REML"))
    
    df <- Cor2DataFrame(x=Cov, n=n, acov=acov)

    ## Effect sizes
    y_wide  <- cbind(df$data[, df$ylabels, drop=FALSE], cluster=cluster)
    
    y_long <- reshape(y_wide, direction = "long", varying=df$ylabels, v.names="y",
                      timevar="r", times=df$ylabels)

    ## FIXME: is it safe?
    y_long <- y_long[order(y_long$id), ]

    ## index <- c(t(outer(seq_len(nrow(df$data)), df$ylabels,
    ##                    function(x,y) paste(x, y, sep="."))))
    
    V_wide <- df$data[, df$vlabels, drop=FALSE]
    V_long <- apply(V_wide, 1, function(x) vec2symMat(unlist(x), diag=TRUE, byrow=FALSE),
                    simplify=FALSE)
    V_long <- bdiagMat(V_long)

    ## Remove rows in y with NA and rows and columns in V with NA
    index.NA <- is.na(y_long$y)   
    y_long <- y_long[!index.NA, , drop=FALSE]
    V_long <- V_long[!index.NA, !index.NA, drop=FALSE]

    switch(RE.type,
           Diag = struct <- "DIAG",
           Symm = struct <- "UN",
           Zero = struct <- "ID")

    fit <- eval(parse(text="metafor::rma.mv(y, V_long, mods = ~ r-1, random = ~ r | id,
                           struct=struct, method=method, sparse=TRUE, data=y_long, ...)"))

    out <- list(rma.fit=fit, cluster=y_long$cluster, n=n, obslabels=rownames(Cov[[1]]), ylabels=df$ylabels)    
    class(out) <- "tssemRobust1"
    out
}

coef.tssemRobust1 <- function(object, ...) {
    if (!is.element("tssemRobust1", class(object)))
        stop("\"object\" must be an object of class \"tssemRobust1\".")
    out <- coef(object$rma.fit, ...)

    ## Arrange the output according to the order of ylabels
    out <-  out[paste0("r", object$ylabels)]
    names(out) <- object$ylabels
    out
}

vcov.tssemRobust1 <- function(object, ...) {
    if (!is.element("tssemRobust1", class(object)))
        stop("\"object\" must be an object of class \"tssemRobust1\".")
    vc <- metafor::robust(object$rma.fit, cluster=object$cluster, ...)
    out <- vc$vb

    ## Arrange the output according to the order of ylabels
    labels <- paste0("r", object$ylabels)
    out <-  out[labels, labels]
    dimnames(out) <- list(object$ylabels, object$ylabels)
    out    
}

summary.tssemRobust1 <- function(object, ...) {
    metafor::robust(object$rma.fit, cluster=object$cluster, ...)
}

tssemRobust2 <- function(tssem1.obj, RAM=NULL, Amatrix=NULL, Smatrix=NULL, Fmatrix=NULL,
                         diag.constraints=FALSE, intervals.type = c("z", "LB"),
                         mx.algebras=NULL, mxModel.Args=NULL,
                         model.name=NULL, suppressWarnings=TRUE, silent=TRUE, run=TRUE,
                         adjust=FALSE, ...) {
    if ( !is.element( class(tssem1.obj), "tssemRobust1") )
        stop("\"tssem1.obj\" must be of class \"tssemRobust1\".")

    Cov <- vec2symMat(coef(tssem1.obj), diag=FALSE)
    dimnames(Cov) <- list(tssem1.obj$obslabels, tssem1.obj$obslabels)

    wls(Cov=Cov, aCov=vcov(tssem1.obj, adjust=adjust), n=sum(tssem1.obj$n),
        RAM=RAM, Amatrix=Amatrix, Smatrix=Smatrix, Fmatrix=Fmatrix,
        diag.constraints=diag.constraints, cor.analysis=TRUE,
        intervals.type=intervals.type, mx.algebras=mx.algebras, mxModel.Args=mxModel.Args,
        model.name=model.name, suppressWarnings=suppressWarnings,
        silent=silent, run=run, ...) 
}

