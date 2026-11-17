#' Summary Estimate Correlation (SEC) of Meta-Analyses with Overlapping
#' Primary Studies
#'
#' It computes the summary estimate correlation (SEC) matrix of several
#' meta-analyses from their citation (incidence) matrix, together with the
#' corrected covered area (CCA). When the variances of the summary estimates
#' are given, it also returns their approximate sampling covariance matrix,
#' which can be used in a second-order meta-analysis.
#'
#' Let \eqn{X} be the \eqn{S \times m}{S x m} citation matrix, where
#' \eqn{X_{si}=1}{X[s,i]=1} when primary study \eqn{s} is included in
#' meta-analysis \eqn{i}. SEC is the correlation between two meta-analytic
#' summary estimators under a specified covariance model. Three versions are
#' available:
#' \describe{
#' \item{\code{"count"}}{The count-based approximation
#' \eqn{R_{ij}=k_{ij}/\sqrt{k_i k_j}}{R[i,j]=k[i,j]/sqrt(k[i]*k[j])}, where
#' \eqn{k_i}{k[i]} is the number of studies in meta-analysis \eqn{i} and
#' \eqn{k_{ij}}{k[i,j]} is the number of studies shared by meta-analyses
#' \eqn{i} and \eqn{j}. It requires only \code{X}. It is exact for unweighted
#' means of independent study estimates with equal variances and approximate
#' otherwise.}
#' \item{\code{"size"}}{The size-weighted approximation
#' \eqn{R_{ij}=\sum_{s \in i \cap j} n_s/\sqrt{\sum_{s \in i} n_s \sum_{s \in
#' j} n_s}}{R[i,j]=sum(n[shared])/sqrt(sum(n[i])*sum(n[j]))}, where
#' \eqn{n_s}{n[s]} is the sample size of primary study \eqn{s}. It assumes
#' that the sampling variances are inversely proportional to the sample sizes
#' and that both meta-analyses use inverse-variance weights. It ignores the
#' heterogeneity variance and therefore approximates the conditional target.}
#' \item{\code{"weighted"}}{The exact weighted SEC for fixed coefficients. Let
#' \eqn{C} be the matrix of the normalized coefficients, with
#' \eqn{C_{si}=X_{si}w_{si}/\sum_t X_{ti}w_{ti}}{C[s,i]=X[s,i]*w[s,i]/sum(X[,i]*w[,i])},
#' and \eqn{Q=\mathrm{Diag}(q_1,\ldots,q_S)}{Q=Diag(q)}. The covariance matrix
#' of the summary estimates is \eqn{C^T Q C}{t(C) Q C}, and \eqn{R} is its
#' correlation matrix. Use the sampling variances \eqn{v_s}{v[s]} as
#' \eqn{q_s}{q[s]} for the conditional target and \eqn{v_s+\tau^2}{v[s]+tau2}
#' for the marginal target. When the weights are estimated from the same data,
#' e.g., with a REML estimate of \eqn{\tau^2}{tau2}, the result is a plug-in
#' approximation.}
#' }
#' Each row of \code{X} must be a distinct and independent study-level unit
#' that contributes the same estimand to every meta-analysis that includes it.
#' A paper cited by two meta-analyses that use different outcomes, time
#' points, subgroups, or contrasts is not a shared study.
#'
#' The CCA and the pairwise CCA are computed with
#' \code{\link[metaSEM]{CCA}}. CCA describes the overlap of the inclusion
#' sets; it is not the correlation between the summary estimates.
#'
#' If \code{V} is given, the approximate sampling covariance matrix of the
#' summary estimates is \eqn{\Sigma=D_V R D_V}{Sigma = D R D}, where
#' \eqn{D_V}{D} is a diagonal matrix of \eqn{\sqrt{V_i}}{sqrt(V)}. \code{V} and
#' \eqn{R} should refer to the same target, e.g., the variances of
#' random-effects summaries with the marginal target. The SEC matrix is
#' singular when the effective study coefficient vectors are linearly
#' dependent; with column-specific weights this need not coincide with
#' dependence of the columns of \code{X}. A warning is given for a singular
#' or nearly singular matrix; a redundant meta-analysis may be removed before
#' the second-order analysis.
#'
#' @param X A numeric or logical matrix or a data frame of 0 and 1, with the
#' primary studies in the rows and the meta-analyses in the columns. Missing
#' values are treated as 0 (not included). Rows without any inclusion are
#' removed. Associated \code{n}, \code{q}, and weight values for removed rows
#' are ignored; values for included rows must be finite.
#' @param method Either \code{"count"} (the default), \code{"size"}, or
#' \code{"weighted"}. See Details.
#' @param n A vector of the sample sizes of the primary studies, one for each
#' row of \code{X}. It is required when \code{method="size"}.
#' @param q A vector of the variances of the primary-study estimates, one for
#' each row of \code{X}. It is required when \code{method="weighted"}.
#' @param weights Unnormalized weights of the primary studies when
#' \code{method="weighted"}. It is either a vector, with one weight for each
#' row of \code{X}, used in all meta-analyses, or a matrix of the same
#' dimensions as \code{X}. If it is missing, the inverse-variance weights
#' \code{1/q} are used.
#' @param V An optional vector of the sampling variances of the summary
#' estimates, one for each column of \code{X}.
#' @return An object of class \code{SEC}, a list with the components:
#' \item{R}{The SEC matrix.}
#' \item{Sigma}{The approximate sampling covariance matrix of the summary
#' estimates, or \code{NULL} if \code{V} is not given.}
#' \item{k}{The numbers of studies in the meta-analyses.}
#' \item{n.studies}{The number of distinct studies in \code{X}.}
#' \item{overlap}{The matrix \eqn{X^T X}{t(X) X} with the numbers of studies
#' on the diagonal and the numbers of shared studies off the diagonal.}
#' \item{CCA}{The corrected covered area of all meta-analyses. See
#' \code{\link[metaSEM]{CCA}}.}
#' \item{CCA.pairwise}{The matrix of the pairwise corrected covered areas.}
#' \item{eigenvalues}{The eigenvalues of \code{R}.}
#' \item{method}{The method used.}
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @references Cheung, M. W.-L. (2026). \emph{Summary estimate correlation
#' (SEC): Quantifying overlap-induced dependence in meta-analyses}.
#' Manuscript in preparation.
#'
#' Pieper, D., Antoine, S.-L., Mathes, T., Neugebauer, E. A. M.,
#' & Eikermann, M. (2014). Systematic review finds overlapping reviews were
#' not mentioned in every other overview. \emph{Journal of Clinical
#' Epidemiology}, \bold{67}(4), 368--375. \doi{10.1016/j.jclinepi.2013.11.007}
#' @seealso \code{\link[metaSEM]{CCA}}, \code{\link[metaSEM]{Turk23}}
#' @keywords meta-analysis
#' @examples
#' ## A citation matrix of 4 primary studies and 3 meta-analyses
#' X <- matrix(c(1, 1, 0,
#'               1, 1, 0,
#'               1, 0, 1,
#'               0, 0, 1), ncol=3, byrow=TRUE,
#'             dimnames=list(paste0("Study", 1:4), c("MA1", "MA2", "MA3")))
#' SEC(X)
#'
#' ## Size-weighted SEC with the sample sizes of the primary studies
#' SEC(X, method="size", n=c(50, 100, 200, 80))
#'
#' ## Weighted SEC with the sampling variances of the primary studies.
#' ## Since q = 4/n, it equals the size-weighted SEC.
#' SEC(X, method="weighted", q=c(0.08, 0.04, 0.02, 0.05))
#'
#' ## Second-order meta-analysis of the data in Turk23
#' ## Standard errors of the summary estimates from their 95% CIs
#' se <- with(Turk23$data, (upper-lower)/(2*qnorm(0.975)))
#' SEC(Turk23$X, V=se^2)
SEC <- function(X, method=c("count", "size", "weighted"), n, q, weights, V) {
  method <- match.arg(method)

  X <- .citation.matrix(X)

  keep <- rowSums(X) > 0
  ## Values for excluded studies do not contribute to the calculation.
  check.vec <- function(x, name) {
    if (!is.numeric(x) || length(x) != nrow(X) || any(!is.finite(x[keep])))
      stop("'", name, "' must be a numeric vector with finite values for included studies, one for each row of 'X'.\n")
    ## Drop names and dimensions, e.g., of a one-dimensional array from tapply()
    as.numeric(x)
  }
  if (!any(keep)) stop("'X' does not include any study.\n")

  if (method=="size") {
    if (missing(n)) stop("'n' is required when method=\"size\".\n")
    n <- check.vec(n, "n")
    if (any(n[keep] <= 0)) stop("'n' must be positive.\n")
  } else if (method=="weighted") {
    if (missing(q)) stop("'q' is required when method=\"weighted\".\n")
    q <- check.vec(q, "q")
    if (any(q[keep] <= 0)) stop("'q' must be positive.\n")
    if (missing(weights)) {
      weights <- numeric(length(q))
      weights[keep] <- 1/q[keep]
    } else if (is.matrix(weights) || is.data.frame(weights)) {
      weights <- as.matrix(weights)
      if (!is.numeric(weights) || !identical(dim(weights), dim(X)) ||
          any(!is.finite(weights[keep,,drop=FALSE])))
        stop("A matrix of 'weights' must have the same dimensions as 'X' and be finite for included studies.\n")
    } else {
      weights <- check.vec(weights, "weights")
    }
    included.weights <- if (is.matrix(weights)) weights[keep,,drop=FALSE] else weights[keep]
    if (any(included.weights < 0)) stop("'weights' must be nonnegative for included studies.\n")
    weights <- matrix(weights, nrow=nrow(X), ncol=ncol(X))
  }

  X <- X[keep, , drop=FALSE]
  k <- colSums(X)
  if (any(k == 0))
    stop("These meta-analyses do not include any study: ",
         paste(colnames(X)[k == 0], collapse=", "), ".\n")

  overlap <- crossprod(X)
  S <- nrow(X)
  m <- ncol(X)

  R <- switch(method,
              count = cov2cor(overlap),
              size = cov2cor(crossprod(X, n[keep]*X)),
              weighted = {
                C <- X*weights[keep, , drop=FALSE]
                Csum <- colSums(C)
                if (any(Csum == 0))
                  stop("The weights of these meta-analyses sum to 0: ",
                       paste(colnames(X)[Csum == 0], collapse=", "), ".\n")
                C <- sweep(C, 2, Csum, "/")
                cov2cor(crossprod(C, q[keep]*C))
              })
  dimnames(R) <- list(colnames(X), colnames(X))

  eigenvalues <- eigen(R, symmetric=TRUE, only.values=TRUE)$values
  if (min(eigenvalues) < 1e-8*max(eigenvalues))
    warning("The SEC matrix is singular or nearly singular; its weighted study contributions may be linearly dependent. ",
            "Consider removing a redundant meta-analysis.\n")

  if (missing(V)) {
    Sigma <- NULL
  } else {
    if (!is.numeric(V) || length(V) != m || any(!is.finite(V)) || any(V <= 0))
      stop("'V' must be a vector of finite positive variances, one for each column of 'X'.\n")
    D <- diag(sqrt(V), nrow=m)
    Sigma <- D %*% R %*% D
    dimnames(Sigma) <- dimnames(R)
  }

  out <- list(R=R, Sigma=Sigma, k=k, n.studies=S, overlap=overlap, CCA=CCA(X),
              CCA.pairwise=CCA(X, pairwise=TRUE), eigenvalues=eigenvalues,
              method=method)
  class(out) <- "SEC"
  out
}

print.SEC <- function(x, digits=4, ...) {
  if (!inherits(x, "SEC"))
    stop("\"x\" must be an object of class \"SEC\".")
  cat("Method:", switch(x$method, count="count-based SEC",
                        size="size-weighted SEC", weighted="weighted SEC"), "\n")
  cat("Number of distinct studies:", x$n.studies, "\n")
  cat("Corrected covered area (CCA):", round(x$CCA, digits), "\n\n")
  cat("Number of studies (diagonal) and shared studies (off-diagonal):\n")
  print(x$overlap, ...)
  cat("\nPairwise CCA:\n")
  print(round(x$CCA.pairwise, digits), ...)
  cat("\nSummary estimate correlation (SEC) matrix:\n")
  print(round(x$R, digits), ...)
  if (!is.null(x$Sigma)) {
    cat("\nSampling covariance matrix of the summary estimates:\n")
    print(signif(x$Sigma, digits), ...)
  }
  if (min(x$eigenvalues) < 1e-8*max(x$eigenvalues))
    cat("\nThe SEC matrix is singular or nearly singular.\n")
  invisible(x)
}
