list2matrix <- function(x, diag=FALSE) {
  if (!is.list(x))
    stop("\"x\" has to be a list.")
  if (!identical(0, var(sapply(x, function(x){dim(x)[[1]]}))))
    stop("Dimensions of matrices in \"x\" have to be the same in order to stack them together.")
  
  if (is.null(dimnames(x[[1]]))) {
    oldNames <- paste("x", 1:dim(x[[1]])[[1]], sep = "")
  } else {
    oldNames <- dimnames(x[[1]])[[1]]
  }
   
  if (diag) {
    psNames <- vech(outer(oldNames, oldNames, paste, sep = "_"))
    ## out <- t(sapply(x, function(x) {(vech(x))}))
    ## out is a list
    out <- lapply(x, vech)    
  } else {
    psNames <- vechs(outer(oldNames, oldNames, paste, sep = "_"))
    ## Fix a bug found by Steffen Zitzmann when x is a 2x2 matrix with diag=FALSE
    ## It returns 1xn vector rather than nx1 because of sapply().
    ## out <- t(sapply(x, function(x) {(vechs(x))}))
    out <- lapply(x, vechs)
  }

  ## convert the list into a matrix
  ## out <- matrix(unlist(out), nrow=length(out), byrow=TRUE)
  out <- do.call(rbind, out)
  
  dimnames(out) <- list(names(x), psNames)
  out
}
