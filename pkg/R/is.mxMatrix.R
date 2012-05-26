is.mxMatrix <- function(x, ...) {
  test <- class(x)
  if (length(grep("Matrix", test[1]))==0) {
    return(FALSE)
  } else {
    if (attr(test, "package")=="OpenMx") {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}
