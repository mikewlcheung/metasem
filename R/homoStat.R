#' Test the Homogeneity of Effect Sizes
#' 
#' It tests the homogeneity of univariate and multivariate effect sizes.
#' 
#' 
#' @param y A vector of effect size for univariate meta-analysis or a
#' \eqn{k}{k} x \eqn{p}{p} matrix of effect sizes for multivariate
#' meta-analysis where \eqn{k}{k} is the number of studies and \eqn{p}{p} is
#' the number of effect sizes.
#' @param v A vector of the sampling variance of the effect size for univariate
#' meta-analysis or a \eqn{k}{k} x \eqn{p*}{p*} matrix of the sampling
#' covariance matrix of the effect sizes for multivariate meta-analysis where
#' \eqn{p* = p(p+1)/2 }{p* = p(p+1)/2}. It is arranged by column major as used
#' by \code{\link[OpenMx]{vech}}. It is assumed that there is no missing value
#' in \code{v} if \code{y} is complete. If there are missing values in \code{v}
#' due to the missingness on \code{y}, the missing values in \code{v} will be
#' removed automatically.
#' @return A list of \item{Q}{Q statistic on the null hypothesis of homogeneity
#' of effect sizes. It has an approximate chi-square distribution under the
#' null hypothesis.} \item{Q.df}{Degrees of freedom of the Q statistic}
#' \item{pval}{p-value on the test of homogeneity of effect sizes}
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{meta}}
#' @references Becker, B. J. (1992). Using results from replicated studies to
#' estimate linear models. \emph{Journal of Educational Statistics}, \bold{17},
#' 341-362.
#' 
#' Cheung, M. W.-L. (2010). Fixed-effects meta-analyses as multiple-group
#' structural equation models. \emph{Structural Equation Modeling}, \bold{17},
#' 481-509.
#' 
#' Cochran, W. G. (1954). The combination of estimates from different
#' experiments. \emph{Biometrics}, \bold{10}, 101-129.
#' @keywords meta-analysis
#' @examples
#' 
#' with( Hox02, homoStat(yi, vi) )
#' 
#' with( HedgesOlkin85, homoStat(y=cbind(d_att, d_ach),
#'       v=cbind(var_att, cov_att_ach, var_ach)) )
#' 
homoStat <- function(y, v) {
  if (is.vector(y)) no.y <- 1 else no.y <- ncol(y)  
  if (is.vector(v)) no.v <- 1 else no.v <- ncol(v)
  if ( no.v != no.y*(no.y+1)/2 )
    stop(paste("The expected no. of columns in v is ", no.y*(no.y+1)/2,
               " while the observed no. of columns in v is ", no.v, ".", sep=""))
    
  if (no.y==1) {
    miss.index <- is.na(y)
    y <- y[!miss.index]
    v <- v[!miss.index]
    w <- 1/v
    beta <- sum(y*w)/sum(w)
    Q <- sum( w*(y-beta)^2 )
    Q.df <- length(y)-1
    pval <- 1-pchisq(Q, df=Q.df)
  } else {
    Y <- matrix( c(t(y)), ncol=1 )
    miss.index <- is.na(Y)
    Y <- matrix( Y[!miss.index], ncol=1 )
    X <- matrix( rep(Diag(no.y), nrow(y)), ncol=no.y, byrow=TRUE )

    X <- X[!miss.index, , drop=FALSE]
    V <- matrix2bdiag(v)
    V <- V[!miss.index, !miss.index, drop=FALSE]
    
    V_inv <- chol2inv(chol(V))
    Q <- t(Y) %*% ( V_inv - V_inv %*% X %*% solve(t(X)
              %*% V_inv %*% X) %*% t(X) %*% V_inv ) %*% Y
    Q.df <- nrow(X)-ncol(X)
	pval <- 1-pchisq(Q, df=Q.df)
  }
    list(Q=Q, Q.df=Q.df, pval=pval)
}


