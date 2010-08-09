\name{asyCov}
\alias{asyCov}
\title{Estimate asymptotic covariance matrix of correlation/covariance matrix
}
\description{It estimates an asymptotic covariance matrix of the correlation/covariance matrix based on multivariate normality. 
}
\usage{
asyCov(x, n, cor.analysis = TRUE, silent = TRUE, suppressWarnings = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A correlation/covariance matrix or a list of correlation/covariance matrices. \code{NA} on the variables or other values defined in \code{na.strings} will be removed before the analysis. Please note that it only checks the diagonal elements of the matrices. If there are missing data, make sure that the diagonals are coded with \code{NA} or other values defined in \code{na.string}.
}
  \item{n}{Sample size or a vector of sample sizes}
  \item{cor.analysis}{Logical. The output is either a correlation or covariance matrix.}
  \item{silent}{Logical. Argument to be passed to \code{\link[OpenMx]{mxRun}}}
  \item{suppressWarnings}{Logical. Argument to be passed to \code{\link[OpenMx]{mxRun}}}
  \item{\dots}{Futher arguments to be passed to \code{\link[OpenMx]{mxRun}}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Asymptotic covariance matrix of the vectorized correlation/covariance matrix}
\references{
Cheung, M. W. L., & Chan, W. (2004). Testing dependent correlation coefficients via structural equation modeling. \emph{Organizational Research Methods}, \bold{7}, 206-223. 
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
C1 <- matrix(c(1,0.5,0.4,0.5,1,0.2,0.4,0.2,1), ncol=3)  
asyCov(C1, n=100)

C2 <- matrix(c(1,0.4,NA,0.4,1,NA,NA,NA,NA), ncol=3)  
C3 <- matrix(c(1,0.2,0.2,1), ncol=2)

asyCov(list(C1,C2,C3), n=c(100,50,50))



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
