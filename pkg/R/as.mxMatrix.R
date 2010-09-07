as.mxMatrix <- function(x, ...) {
  if (!is.matrix(x))
    stop("\"x\" must be a matrix.")
  # suppress warnings
  warn <- options()$warn
  options(warn=-1)
  nRow <- nrow(x)
  nCol <- ncol(x)
  values <- as.numeric(x)  # They are NA for characters 
  free <- is.na(values)    # They are TRUE for parameters with labels 
  freePara1 <- x[free]     # Extract free parameters
  # check if there are any free parameters
  if (length(freePara1)>0) {
    freePara2 <- strsplit(freePara1, "*", fixed=TRUE)
    # replace NA with starting values 0.5 before "0.5*a"
    values[free] <- sapply(freePara2, function(x){ as.numeric(x[1])})
    labels <- matrix(NA, ncol=nCol, nrow=nRow)
    labels[free] <- sapply(freePara2, function(x){ x[2]})
    out <- mxMatrix(type = "Full", nrow=nRow, ncol=nCol, values=values, free=free,
                    labels=labels, ...)
  } else {
    out <- mxMatrix(type = "Full", nrow=nRow, ncol=nCol, values=values, free=free, ...)    
  }
  options(warn=warn)
  out
}
