as.mxMatrix <- function(x, name, ...) {
    ## If it is a vector, the output is a column matrix.
    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }

    # suppress warnings
    ## warn <- options()$warn
    ## options(warn=-1)
    nRow <- nrow(x)
    nCol <- ncol(x)

    # check if "name" was give
    # if not, use the matrix name
    if (missing(name)) {
        name <- as.character(match.call())[2]

        ## Check if "$" is present
        ## Suppose RAM$F, the output name is "F"
        if (grepl("$", name, fixed=TRUE)) {
            name <- strsplit(name, "$", fixed=TRUE)[[1]][2]
        }    
    }

    values <- suppressWarnings(as.numeric(x))  # They are NA for characters
    free <- is.na(values)    # They are TRUE for parameters with labels
    freePara1 <- x[free]     # Extract free parameters
    # check if there are any free parameters

    if (length(freePara1)>0) {
        freePara2 <- strsplit(freePara1, "*", fixed=TRUE)
        # replace NA with starting values 0.5 before "0.5*a"
        values[free] <- suppressWarnings(sapply(freePara2, function(x){ as.numeric(x[1])}))
        labels <- matrix(NA, ncol=nCol, nrow=nRow)
        labels[free] <- sapply(freePara2, function(x){ x[2]})
    
        ## Replace TRUE by FALSE in free when there are definition variables or [,]
        free[grep("data.", labels)] <- FALSE
        free[grep("\\[|,|\\]", labels)] <- FALSE

        ## Check any "@"
        x_pos <- grep("@", x, fixed=TRUE)
        if (length(x_pos)>0) {
            x_at <- strsplit(x=x[x_pos], split="@", fixed=TRUE)
            for (i in seq_along(x_at)) {
                values[x_pos[i]] <- as.numeric(x_at[[i]][1])
                labels[x_pos[i]] <- x_at[[i]][2]
                free[x_pos[i]] <- FALSE                
            }
        }
        
        out <- mxMatrix(type = "Full", nrow=nRow, ncol=nCol, values=values, free=free,
                        name=name, labels=labels, ...)
    } else {
        out <- mxMatrix(type = "Full", nrow=nRow, ncol=nCol, values=values, free=free,
                        name=name, ...)
    }

    ## Add the dimnames only when there are dimnames
    if (!is.null(dimnames(x))) {
        if (!is.null(rownames(x))) {
            if (rownames(x)[1] != "1") {
                ## Make the names valid for the Mmatrix, which has "1" as the rownames
                dim.names <- lapply(dimnames(x), make.names)
                dimnames(out@values) <- dimnames(out@labels) <- dimnames(out@free) <- dim.names
            }
        }
   }
    
  ## options(warn=warn)
  out
}

as.symMatrix <- function(x) {
    if (is.list(x)) {
        for (i in seq_along(x)) {
            x[[i]][] <- vapply(x[[i]], function(z) gsub(".*\\*", "", z), character(1))
        }
    } else {
        x[] <- vapply(x, function(z) gsub(".*\\*", "", z), character(1))
    }
    x
}

as.mxAlgebra <- function(x, startvalues=NULL, name="X") {
    ## If it is a vector, the output is a column matrix.
    if (!is.matrix(x)) {
        x <- as.matrix(x)
    }
    
    ## Xvars: a column vector of the parameters in x, e.g., a, b, c.
    vars <- sort(all.vars(parse(text=x)))
    ## Provide 0 as starting value 
    Xvars <- create.mxMatrix(paste0("0*", vars), ncol=1, nrow=length(vars))
    ## Change the matrix name
    Xvars@name <- paste0(name, "vars")

    ## Change the starting values in Xvars, if provided
    if (!is.null(startvalues)) {
        for (i in seq_along(startvalues)) {
            starti <- unlist(startvalues[i])
            Xvars$values[Xvars$labels==names(starti)] <- starti
        }
    }
    
    xrow <- nrow(x)
    xcol <- ncol(x)

    ## This is general but not efficient
    ## The mxAlgebra matrix is built on rbind(cbind(...)) of each Xi_j
    ## If x[i,j] = a^b, then Xij <- mxAlgebra(a^b, name="xi_j")
    for (j in seq_len(xcol))
        for (i in seq_len(xrow)) {
            ## Elements of each Xi_j
            Xij <- paste0(name,i,"_",j, " <- mxAlgebra(",x[i,j], ", name='",name,i,"_",j,"')")
            eval(parse(text=Xij))               
    }

    ## Xmat: a matrix of the named Xi_j of mxalgebra
    Xmat <- outer(seq_len(xrow), seq_len(xcol), function(x, y) paste0(name,x,"_",y))
    Xmat <- paste0("cbind(", apply(Xmat, 1, paste0, collapse=", "), ")")
    Xmatrix  <- paste0(name, " <- mxAlgebra(rbind(", paste0(Xmat, collapse=", "),
                       "), name='", name,"')")
    eval(parse(text=Xmatrix))

    ## Prepare mxalgebra matrices for output
    ## "name", e.g., A: the mxAlgebra object
    ## Xvars, e.g., Avars: the variables (parameters) to build the mxAlgebra
    ## names: a list of the names of these matrices
    ## Xlist, e.g., Alist: a list of mxalgebra of X1_1, X1_2, X2_1, ...
    Xlist <- outer(seq_len(xrow), seq_len(xcol), function(x,y) paste0(name,x,"_",y))
    Xnames <- c(name, Xvars@name, Xlist)
    Xlist <- paste0("list(", paste0(Xlist, "=", Xlist, collapse=", "), ")")
    Xlist <- eval(parse(text=Xlist))
    # list(Xmatrix=Xmatrix, Xvars=Xvars, Xnames=Xnames, Xlist=Xlist)

    ## out <- paste0("out <- list(", name, "=", name, ", names=Xnames, ", name, 
    ##               "vars=Xvars, ", name, "list=Xlist",")" )
    out <- paste0("out <- list(mxalgebra=", name, ", parameters=Xvars, list=Xlist",")" )
    eval(parse(text=out))
    out
}
