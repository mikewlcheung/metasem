checkRAM <- function(Amatrix, Smatrix, cor.analysis=TRUE) {
    if (missing(Amatrix)&missing(Smatrix)) {
        warning("Either 'Amatrix' or 'Smatrix' must be present.")
    }    

    ## Check A
    if (!missing(Amatrix)) {
        ## Convert A into mxMatrix if they are not yet.
        if (is.matrix(Amatrix)) {
            Amatrix <- as.mxMatrix(Amatrix, name="A")
        } else {
            Amatrix@name <- "A"
        }

        ## A_name <- deparse(substitute(A))

        ## Check diagonals: either free or non-zero values
        if ( any(Diag(Amatrix$free)) | any(Diag(Amatrix$values)!=0) ) {
            warning("Diagonals of the 'Amatrix' must be zeros.\n")
        }

        ## Check both lower tri & upper tri are TRUE
        if ( any(Amatrix$free & t(Amatrix$free)) ) {
            warning("Non-recursive models are not allowed in the 'Amatrix'.\n")
        }
    }

    ## Check S
    if (!missing(Smatrix)) {
        if (is.matrix(Smatrix)) {
            Smatrix <- as.mxMatrix(Smatrix, name="S")
        } else {
            Smatrix@name <- "S"
        }

        ## S_name <- deparse(substitute(S))
        ## Cannot check 'free' for definition variables
        ## Only check labels!!!
        S_labels <- Smatrix$labels

        ## Check symmetric
        if (!isSymmetric(Smatrix$free)) {
            warning("The free parameters of the 'Smatrix' must be symmetric.\n")
        }

        if (!isSymmetric(Smatrix$labels)) {
            warning("The labels of 'Smatrix' must be symmetric.\n")
        }

        if (!isSymmetric(Smatrix$values)) {
            warning("The values of 'Smatrix' must be symmetric.\n")
        }
        
        ## ## Check digaonals
        ## if (SdiagZero==TRUE & any(!is.na(Diag(S_labels)))) {
        ##     warning("Diagonal of 'S' must be 0.\n")
        ## }

        ## Check both A and S: Variances of IVs must be fixed at 1 and DVs must be free
        ## Limitation: it may still give warnings when there are DVs with fixed parameters in A
        if (cor.analysis==TRUE & !missing(Amatrix)) {
            ## Check A if it is a DV
            dv <- apply(Amatrix$free, 1, any)

            if ( any(diag(Smatrix$free)[!dv]) | !all(diag(Smatrix$values)[!dv]==1) ) {
                warning("The variances of the independent variables in 'Smatrix' must be fixed at 1.")
            }            
            if (!all(diag(Smatrix$free)[dv])) {
                warning("The variances of the dependent variables in 'Smatrix' should be free.")
            }                        
        }        
    }      
}        
