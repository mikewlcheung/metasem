\name{list2matrix}
\alias{list2matrix}
\title{Convert a List of Symmetric Matrices into a Stacked Matrix
}
\description{It converts a list of symmetric matrices into
  a stacked matrix. Dimensions of the symmetric matrices have to be the
  same. It tries to preserve the dimension names if possible. Dimension names will
  be created if there is no dimension names in the first symmetric matrix. 
}
\usage{
list2matrix(x, diag = FALSE)
}
\arguments{
  \item{x}{A list of \eqn{k}{k} \eqn{p}{p} x \eqn{p}{p} symmetric matrices.
}
  \item{diag}{Logical. If it is \code{TRUE}, \code{\link[OpenMx]{vech}}
	is used to vectorize the (covariance) matrices. If it is \code{FALSE}, \code{\link[OpenMx]{vechs}}
	is used to vectorize the (correlation) matrices.
}
}
\details{
}
\value{A \eqn{k}{k} x \eqn{p*}{p*} stacked matrix where \eqn{p* = p(p-1)/2 }{p*
	= p(p-1)/2} for \code{diag}=\code{FALSE} or \eqn{p* = p(p+1)/2
  }{p* = p(p+1)/2} for \code{diag}=\code{TRUE}.

}
\references{
%% ~put references to the literature/web site here ~
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
C2 <- matrix(c(1,0.4,NA,0.4,1,NA,NA,NA,NA), ncol=3)  

## A list without dimension names 
list2matrix(list(C1, C2))
#      x2_x1 x3_x1 x3_x2
# [1,]   0.5   0.4   0.2
# [2,]   0.4    NA    NA

dimnames(C1) <- list( c("x","y","z"), c("x","y","z") )
dimnames(C2) <- list( c("x","y","z"), c("x","y","z") )

## A list with dimension names
list2matrix(list(C1, C2))
#      y_x z_x z_y
# [1,] 0.5 0.4 0.2
# [2,] 0.4  NA  NA
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
