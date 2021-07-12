## Easy creation of mx models
## Note. Dimension names of variables are assumed stored in RAM$F;
## otherwise, they have to be specified in var.names. 
create.mxModel <- function(model.name="mxModel", RAM=NULL,
                           data, intervals.type = c("z", "LB"),
                           var.names=NULL, mxModel.Args=NULL,
                           run=TRUE, mxTryHard=FALSE,
                           silent=TRUE, ...) {

    intervals.type <- match.arg(intervals.type)
    switch(intervals.type,
           z = intervals <- FALSE,
           LB = intervals <- TRUE)

    Amatrix <- as.mxMatrix(RAM$A, name="Amatrix")
    Smatrix <- as.mxMatrix(RAM$S, name="Smatrix")
    Fmatrix <- as.mxMatrix(RAM$F, name="Fmatrix")
    Mmatrix <- as.mxMatrix(RAM$M, name="Mmatrix")
    
    ## Some basic checking in RAM
    checkRAM(Amatrix, Smatrix, cor.analysis=FALSE)

    ## Extract the dimnames from Fmatrix$values
    if (is.null(var.names)) {
        var.names <- colnames(Fmatrix$values)
    }

    ## If matrix or data.frame is provided, setup mxData
    ## Otherwise, users have to setup mxData()
    if (is.data.frame(data) | is.matrix(data)) {
        mx.data <- mxData(observed=data, type="raw")
    } else {
        mx.data <- data
    }    

    ## Create an identity matrix from the dimensions of Amatrix
    Id <- as.mxMatrix(diag(nrow(Amatrix$values)), name="Id")

    ## Expected covariance matrix and means of the observed and latent variables
    Id_A <- mxAlgebra(solve(Id - Amatrix), name="Id_A")
    expCov <- mxAlgebra(Id_A %&% Smatrix, name="expCov")
    expMean <- mxAlgebra(Mmatrix %*% t(Id_A), name="expMean")
    
    mx.model <- mxModel(model.name, Amatrix, Smatrix, Fmatrix, Mmatrix,
                        Id, Id_A, expCov, expMean,
                        mx.data, mxFitFunctionML(),
                        mxCI(c("Amatrix", "Smatrix", "Mmatrix")),
                        mxExpectationRAM(A="Amatrix", S="Smatrix",
                                         F="Fmatrix", M="Mmatrix",
                                         dimnames=var.names))
    
    ## Add additional arguments to mxModel
    if (!is.null(mxModel.Args)) {
        for (i in seq_along(mxModel.Args)) {
            mx.model <- mxModel(mx.model, mxModel.Args[[i]])
        }
    }
    
    ## Add mxAlgebra and mxConstraint from RAM$mxalgebra
    if (!is.null(RAM$mxalgebras)) {
            for (i in seq_along(RAM$mxalgebras)) {
                mx.model <- mxModel(mx.model, RAM$mxalgebras[[i]])
                ## Name of the mxalgebra
                name.mxalgebra <- names(RAM$mxalgebras)[i]
                ## Check if the name constains constraint1, constraint2, ...,
                ## If no, they are mxalgebra, not mxconstraints. Include them in mxCI. 
                if (!grepl("^constraint[0-9]", name.mxalgebra)) {
                    mx.model <- mxModel(mx.model, mxCI(c(name.mxalgebra)))
                }
            }
    }
    
    if (run==TRUE) {
        if (mxTryHard==TRUE) {
            out <- mxTryHard(mx.model, intervals=intervals, silent=silent, ...)
        } else {
            out <- mxRun(mx.model, intervals=intervals, silent=silent, ...)
        }
    } else {
        out <- mx.model
    }
    out
}
