## OpenMx uses a period to qualify references. Reject periods in model and
## data names before they can be mistaken for such references.
.check.no.period <- function(x, what) {
  x <- unique(as.character(x[!is.na(x)]))
  bad <- x[grepl(".", x, fixed=TRUE)]
  if (length(bad)) {
    example <- paste0("'", bad[1], "' to '",
                      gsub(".", "_", bad[1], fixed=TRUE), "'")
    stop(what, " contain unsupported names with periods: ",
         paste(sprintf("'%s'", bad), collapse=", "),
         ". Rename consistently in the model, data, matrix dimnames and ",
         "subset.variables, for example ", example,
         "; check that the new names are unique.", call.=FALSE)
  }
  invisible(NULL)
}

.check.parameter.names <- function(labels, free=NULL) {
  labels <- as.character(labels)
  if (is.null(free)) free <- rep(TRUE, length(labels))
  definition <- !is.na(labels) & startsWith(labels, "data.")
  .check.no.period(substring(labels[definition], 6L), "Definition-variable column names")
  ## Fixed qualified references such as submodel.parameter[1,1] belong to
  ## OpenMx, rather than being ordinary free parameter labels.
  qualified <- !free & grepl("^[A-Za-z][A-Za-z0-9_]*\\.[A-Za-z][A-Za-z0-9_]*(\\[.*\\])?$", labels)
  .check.no.period(labels[!definition & !qualified], "Parameter labels")
}

.check.model.matrix <- function(x) {
  if (is.null(x)) return(invisible(NULL))
  if (inherits(x, "MxMatrix")) {
    .check.no.period(unlist(dimnames(x$values), use.names=FALSE), "Model variable names")
    .check.parameter.names(x$labels, x$free)
  } else if (is.matrix(x)) {
    .check.no.period(unlist(dimnames(x), use.names=FALSE), "Model variable names")
    .check.matrix.labels(x)
  }
  invisible(NULL)
}

.check.matrix.labels <- function(x) {
  if (is.character(x)) {
    entries <- x[!is.na(x)]
    entries <- entries[is.na(suppressWarnings(as.numeric(entries)))]
    fixed <- grepl("@|\\[", entries)
    labels <- sub("^.*[*@]", "", entries)
    .check.parameter.names(labels, !fixed)
  }
  invisible(NULL)
}

.check.ram.names <- function(x) {
  if (is.null(x) || !is.list(x)) return(invisible(NULL))
  if (all(c("A", "S", "F") %in% names(x))) {
    for (part in c("A", "S", "F", "M")) .check.model.matrix(x[[part]])
  } else {
    for (part in x) if (is.list(part)) .check.ram.names(part)
  }
  invisible(NULL)
}

## Recurse over groups/studies/clusters without treating their list names as
## variable names: group and study identifiers may themselves contain dots.
.check.cov.names <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    .check.no.period(unlist(dimnames(x), use.names=FALSE), "Covariance variable names")
  } else if (is.list(x)) {
    for (part in x) .check.cov.names(part)
  }
  invisible(NULL)
}

.check.mean.names <- function(x) {
  if (is.list(x)) {
    for (part in x) .check.mean.names(part)
  } else if (is.matrix(x)) {
    .check.no.period(colnames(x), "Mean variable names")
  } else {
    .check.no.period(names(x), "Mean variable names")
  }
  invisible(NULL)
}

.check.meta.input.names <- function(call, data=NULL) {
  if (!is.null(data)) {
    columns <- if (is.list(data) && !is.data.frame(data)) names(data) else colnames(data)
    .check.no.period(columns, "Data column names")
  }
  ## Names of R objects used in expressions are not model variable names.
  ## Actual data columns and evaluated matrix column names are checked by
  ## the caller; rejecting all.vars() here also rejects e.g. my.matrix[, 1].
}

## Inspect lavaan models used by the wrappers without imposing lavaan2RAM's
## narrower operator or automatic-option restrictions.
.check.lavaan.model.names <- function(model) {
  pt <- tryCatch(suppressWarnings(lavaan::lavaanify(model)), error=function(e) NULL)
  if (!is.null(pt)) {
    paths <- pt$op %in% c("=~", "<~", "~", "~~", "~1", "|")
    .check.no.period(c(pt$lhs[paths], pt$rhs[paths]), "Model variable names")
    .check.parameter.names(pt$label[paths | pt$op == ":="])
  }
  invisible(NULL)
}

## Unknown sampling covariances between retained effects cannot be replaced
## by a large number: that makes the sampling covariance non-positive-definite.
.check.missing.sampling.cov <- function(y, v) {
  if (is.null(dim(y)) || ncol(y) < 2L) return(invisible(NULL))
  p <- ncol(y)
  if (is.null(dim(v)) || ncol(v) != p*(p+1L)/2L) return(invisible(NULL))
  index <- matrix(NA_integer_, p, p)
  index[lower.tri(index, diag=TRUE)] <- seq_len(ncol(v))
  index[upper.tri(index)] <- t(index)[upper.tri(index)]
  for (row in seq_len(nrow(y))) {
    observed <- !is.na(y[row, ]) & !is.na(v[row, diag(index)])
    for (i in seq_len(p)) for (j in seq_len(i-1L))
      if (observed[i] && observed[j] && is.na(v[row, index[i,j]]))
        stop("Missing sampling covariance between observed effects ", j, " and ", i,
             " in study row ", row, ". Supply it or explicitly choose an appropriate approximation.", call.=FALSE)
  }
  invisible(NULL)
}
