#' Corrected Covered Area (CCA) of Meta-Analyses with Overlapping Primary
#' Studies
#'
#' It computes the corrected covered area (CCA) of several meta-analyses, or
#' the pairwise CCA of each pair of meta-analyses, from their citation
#' (incidence) matrix.
#'
#' Let \eqn{X} be the \eqn{S \times m}{S x m} citation matrix, where
#' \eqn{X_{si}=1}{X[s,i]=1} when primary study \eqn{s} is included in
#' meta-analysis \eqn{i}. The CCA (Pieper et al., 2014) is
#' \deqn{\mathrm{CCA}=\frac{N-S}{S(m-1)},}{CCA = (N-S)/(S*(m-1)),}
#' where \eqn{N} is the total number of inclusions, i.e., the sum of
#' \eqn{X}, and \eqn{S} is the number of distinct studies. The numerator
#' counts the inclusions of the studies after their first one. CCA ranges from
#' 0, when no study is included in more than one meta-analysis, to 1, when all
#' studies are included in all meta-analyses. Pieper et al. (2014) suggested
#' that CCA values of 0 to 0.05, 0.06 to 0.10, 0.11 to 0.15, and above 0.15
#' indicate slight, moderate, high, and very high overlap, respectively.
#'
#' The pairwise CCA of meta-analyses \eqn{i} and \eqn{j} is the CCA of these
#' two columns of \eqn{X}, which equals
#' \eqn{k_{ij}/(k_i+k_j-k_{ij})}{k[i,j]/(k[i]+k[j]-k[i,j])}, where
#' \eqn{k_i}{k[i]} is the number of studies in meta-analysis \eqn{i} and
#' \eqn{k_{ij}}{k[i,j]} is the number of studies shared by meta-analyses
#' \eqn{i} and \eqn{j}.
#'
#' CCA describes the overlap of the inclusion sets. It is not the correlation
#' between the summary estimates of the meta-analyses, and the same CCA may
#' come from a few pairs with high overlap or from many pairs with modest
#' overlap. Use \code{\link[metaSEM]{SEC}} for the correlation matrix and the
#' sampling covariance matrix of the summary estimates.
#'
#' @param X A numeric or logical matrix or a data frame of 0 and 1, with the
#' primary studies in the rows and the meta-analyses in the columns. Missing
#' values are treated as 0 (not included). Rows without any inclusion are
#' removed.
#' @param pairwise Logical. If it is \code{FALSE} (the default), the CCA of
#' all meta-analyses is returned. If it is \code{TRUE}, the matrix of the
#' pairwise CCA is returned.
#' @return The CCA of all meta-analyses if \code{pairwise=FALSE}, or a
#' symmetric matrix of the pairwise CCA, with 1 on the diagonal, if
#' \code{pairwise=TRUE}.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @references Pieper, D., Antoine, S.-L., Mathes, T., Neugebauer, E. A. M.,
#' & Eikermann, M. (2014). Systematic review finds overlapping reviews were
#' not mentioned in every other overview. \emph{Journal of Clinical
#' Epidemiology}, \bold{67}(4), 368--375. \doi{10.1016/j.jclinepi.2013.11.007}
#' @seealso \code{\link[metaSEM]{SEC}}, \code{\link[metaSEM]{Turk23}}
#' @keywords meta-analysis
#' @examples
#' ## A citation matrix of 4 primary studies and 3 meta-analyses
#' X <- matrix(c(1, 1, 0,
#'               1, 1, 0,
#'               1, 0, 1,
#'               0, 0, 1), ncol=3, byrow=TRUE,
#'             dimnames=list(paste0("Study", 1:4), c("MA1", "MA2", "MA3")))
#' CCA(X)
#' CCA(X, pairwise=TRUE)
#'
#' ## CCA of 0.125 reported by Turk et al. (2023)
#' CCA(Turk23$X)
#' round(CCA(Turk23$X, pairwise=TRUE), 3)
CCA <- function(X, pairwise=FALSE) {
  X <- .citation.matrix(X)
  X <- X[rowSums(X) > 0, , drop=FALSE]
  if (nrow(X) == 0) stop("'X' does not include any study.\n")
  k <- colSums(X)
  if (any(k == 0))
    stop("These meta-analyses do not include any study: ",
         paste(colnames(X)[k == 0], collapse=", "), ".\n")

  if (pairwise) {
    overlap <- crossprod(X)
    overlap/(outer(k, k, "+") - overlap)
  } else {
    S <- nrow(X)
    (sum(X) - S)/(S*(ncol(X) - 1))
  }
}

## Check a citation matrix: a numeric matrix of 0 and 1 with named columns
.citation.matrix <- function(X) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X) || !(is.numeric(X) || is.logical(X)))
    stop("'X' must be a numeric or logical matrix or a data frame.\n")
  X[is.na(X)] <- 0
  if (!all(X %in% c(0, 1)))
    stop("'X' must contain only 0 and 1 (or NA for not included).\n")
  X <- X + 0
  if (ncol(X) < 2) stop("'X' must have at least two columns (meta-analyses).\n")
  if (is.null(colnames(X))) colnames(X) <- paste0("MA", seq_len(ncol(X)))
  X
}
