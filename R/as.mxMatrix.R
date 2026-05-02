#' Convert a Matrix into MxMatrix-class
#' 
#' It converts a matrix into \code{MxMatrix-class} via \code{mxMatrix}.
#' 
#' If there are non-numeric values in \code{x}, they are treated as the labels
#' of the parameters. If a "*" is present, the numeric value on the left-hand
#' side will be treated as the starting value for a free parameter. If an "@"
#' is present, the numeric value on the left-hand side will be considered as
#' the value for a fixed parameter. If it is a matrix of numeric values, there
#' are no free parameters in the output matrix.
#' 
#' @param x A character or numeric matrix. If \code{x} is not a matrix,
#' \code{as.matrix(x)} is applied first.
#' @param name An optional character string as the name of the MxMatrix object
#' created by mxModel function. If the \code{name} is missing, the name of
#' \code{x} will be used.
#' @param \dots Further arguments to be passed to
#' \code{\link[OpenMx]{mxMatrix}}. It should be noted that \code{type},
#' \code{nrow}, \code{ncol}, \code{values}, \code{free}, \code{name} and
#' \code{labels} will be created automatically. Thus, these arguments except
#' labels should be avoided in \dots
#' @return A \code{\link[OpenMx]{MxMatrix-class}} object with the same
#' dimensions as \code{x}
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[OpenMx]{mxMatrix}},
#' \code{\link[metaSEM]{create.mxMatrix}},
#' \code{\link[metaSEM]{create.Fmatrix}}, \code{\link[metaSEM]{checkRAM}},
#' \code{\link[metaSEM]{lavaan2RAM}}, \code{\link[metaSEM]{as.symMatrix}}
#' @keywords utilities
#' @examples
#' 
#' ## a and b are free parameters with starting values and labels
#' (a1 <- matrix(c(1:4, "5*a", 6, "7*b", 8, 9), ncol=3, nrow=3))
#' #      [,1] [,2]  [,3] 
#' # [1,] "1"  "4"   "7*b"
#' # [2,] "2"  "5*a" "8"  
#' # [3,] "3"  "6"   "9"  
#' 
#' a1 <- as.mxMatrix(a1)
#' 
#' ## a and b are fixed parameters without any labels, name="new2"
#' (a2 <- matrix(1:9, ncol=3, nrow=3))
#' #      [,1] [,2] [,3]
#' # [1,]    1    4    7
#' # [2,]    2    5    8
#' # [3,]    3    6    9
#' 
#' new2 <- as.mxMatrix(a2, name="new2")
#' 
#' ## Free parameters without starting values
#' (a3 <- matrix(c(1:4, "*a", 6, "*b", 8, 9), ncol=3, nrow=3))
#' #      [,1] [,2] [,3]
#' # [1,] "1"  "4"  "*b"
#' # [2,] "2"  "*a" "8" 
#' # [3,] "3"  "6"  "9" 
#' 
#' a3 <- as.mxMatrix(a3, lbound=0)
#' 
#' ## A free parameter without label
#' (a4 <- matrix(c(1:4, "5*", 6, "7*b", 8, 9), ncol=3, nrow=3))
#' #      [,1] [,2] [,3] 
#' # [1,] "1"  "4"  "7*b"
#' # [2,] "2"  "5*" "8"  
#' # [3,] "3"  "6"  "9"  
#' 
#' a4 <- as.mxMatrix(a4)
#' 
#' ## Convert a scalar into mxMatrix object
#' ## "name" is required as "3*a" is not a valid name.
#' (a5 <- as.mxMatrix("3*a", name="a5"))
#' 
#' ## Free and fixed parameters
#' (a6 <- matrix(c(1, "2*a", "3@b", 4), ncol=2, nrow=2))
#' 
#' as.mxMatrix(a6)
#' 
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



#' Convert a Character Matrix with Starting Values to a Character Matrix
#' without Starting Values
#' 
#' It converts a character matrix with starting values to a character matrix
#' without the starting values.
#' 
#' If there are non-numeric values in \code{x}, they are treated as the labels
#' of the free parameters. If a "*" is present, the numeric value on the
#' left-hand side will be treated as the starting value for a free parameter or
#' a fixed value for a fixed parameter. If it is a matrix of numeric values,
#' there are no free parameters in the output matrix. This function removes the
#' starting values and "*" in the matrices.
#' 
#' @param x A character or numeric matrix or a list of character or numeric
#' matrices.
#' @return A character matrix.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{as.mxMatrix}}
#' @keywords utilities
#' @examples
#' 
#' ## a and b are free parameters with starting values and labels
#' (a1 <- matrix(c(1:4, "5*a", 6, "7*b", 8, 9), ncol=3, nrow=3))
#' #      [,1] [,2]  [,3] 
#' # [1,] "1"  "4"   "7*b"
#' # [2,] "2"  "5*a" "8"  
#' # [3,] "3"  "6"   "9"  
#' 
#' (as.symMatrix(a1))
#' #      [,1] [,2] [,3]
#' # [1,] "1"  "4"  "b" 
#' # [2,] "2"  "a"  "8" 
#' # [3,] "3"  "6"  "9"
#' 
as.symMatrix <- function(x) {
  if (is.list(x)) {
    ## for (i in seq_along(x)) {
    ## Exclude mxalgebras, which creates troubles
    for (i in c("A", "S", "F", "M")) {
      x[[i]][] <- vapply(x[[i]], function(z) gsub(".*\\*", "", z), character(1))
    }
  } else {
    x[] <- vapply(x, function(z) gsub(".*\\*", "", z), character(1))
  }
  x
}



#' Convert a Character Matrix into MxAlgebra-class
#' 
#' It converts a character matrix into \code{MxAlgebra} object.
#' 
#' Suppose the name argument is "X", the output is a list of the following
#' elements.
#' 
#' @param x A character or numeric matrix, which consists of valid operators in
#' \code{mxAlgebra}.
#' @param startvalues A list of starting values of the free parameters. If it
#' is not provided, all free parameters are assumed 0.
#' @param lbound A list of lower bound of the free parameters. If it is not
#' provided, all free parameters are assumed \code{NA}.
#' @param ubound A list of upper bound of the free parameters. If it is not
#' provided, all free parameters are assumed \code{NA}.
#' @param name A character string of the names of the objects based on.
#' @return \item{mxalgebra}{An \code{mxAlgebra} object.} \item{parameters}{A
#' column vector \code{mxMatrix} of the free parameters.} \item{list}{A list of
#' mxMatrix to form the \code{mxAlgebra} object.}
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{as.mxMatrix}}, \code{\link[OpenMx]{mxAlgebra}}
#' @keywords utilities
#' @examples
#' 
#' ## a, b, and c are free parameters
#' (A1 <- matrix(c(1, "a*b", "a^b", "exp(c)"), ncol=2, nrow=2))
#' ##      [,1]  [,2]    
#' ## [1,] "1"   "a^b"   
#' ## [2,] "a*b" "exp(c)"
#' 
#' A <- as.mxAlgebra(A1, startvalues=list(a=1, b=2),
#'                   lbound=list(a=0), ubound=list(b=1, c=2),
#'                   name="A")
#' 
#' ## An object of mxAlgebra
#' A$mxalgebra
#' ## mxAlgebra 'A' 
#' ## $formula:  rbind(cbind(A1_1, A1_2), cbind(A2_1, A2_2)) 
#' ## $result: (not yet computed) <0 x 0 matrix>
#' ## dimnames: NULL
#' 
#' ## A matrix of parameters
#' A$parameters
#' ## FullMatrix 'Avars' 
#' 
#' ## $labels
#' ##      [,1]
#' ## [1,] "a" 
#' ## [2,] "b" 
#' ## [3,] "c" 
#' 
#' ## $values
#' ##      [,1]
#' ## [1,] 1   
#' ## [2,] 2   
#' ## [3,] 0   
#' 
#' ## $free
#' ##      [,1]
#' ## [1,] TRUE
#' ## [2,] TRUE
#' ## [3,] TRUE
#' 
#' ## $lbound
#' ##      [,1]
#' ## [1,]    0
#' ## [2,]   NA
#' ## [3,]   NA
#' 
#' ## $ubound
#' ##      [,1]
#' ## [1,]   NA
#' ## [2,]    1
#' ## [3,]    2
#' 
#' ## A list of matrices of elements for the mxAlgebra
#' A$list
#' ## $A1_1
#' ## mxAlgebra 'A1_1' 
#' ## $formula:  1 
#' ## $result: (not yet computed) <0 x 0 matrix>
#' ## dimnames: NULL
#' 
#' ## $A2_1
#' ## mxAlgebra 'A2_1' 
#' ## $formula:  a * b 
#' ## $result: (not yet computed) <0 x 0 matrix>
#' ## dimnames: NULL
#' 
#' ## $A1_2
#' ## mxAlgebra 'A1_2' 
#' ## $formula:  a^b 
#' ## $result: (not yet computed) <0 x 0 matrix>
#' ## dimnames: NULL
#' 
#' ## $A2_2
#' ## mxAlgebra 'A2_2' 
#' ## $formula:  exp(c) 
#' ## $result: (not yet computed) <0 x 0 matrix>
#' ## dimnames: NULL
#' 
as.mxAlgebra <- function(x, startvalues=NULL, lbound=NULL, ubound=NULL,
                         name="X") {
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

  ## Add lbound
  if (!is.null(lbound)) {
    for (i in seq_along(lbound)) {
      index <- Xvars$labels==names(lbound)[i]
      ## NAs are treated as FALSE
      index[is.na(index)] <- FALSE
      Xvars$lbound[index] <- lbound[[i]]
    }
  }

  ## Add ubound
  if (!is.null(ubound)) {
    for (i in seq_along(ubound)) {
      index <- Xvars$labels==names(ubound)[i]
      index[is.na(index)] <- FALSE
      Xvars$ubound[index] <- ubound[[i]]
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
