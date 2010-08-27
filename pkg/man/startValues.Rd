\name{startValues}
\alias{startValues}
\title{Starting values for tssem1()
}
\description{It calculates starting values for \code{\link{tssem1}} with
  unweighted mean on the correlation/covariance matrices. 
}
\usage{
startValues(x, cor.analysis = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A list of correlation or covariance matrices
}
  \item{cor.analysis}{logical. The output is either a correlation matrix
	or a covariance matrix
}
}
\details{It calculates the starting values with the unweighted
  means on the correlation/covariance matrices. Missing values are removed before the calculations.
  \link[Matrix]{nearPD} is used to ensure that the calculated matrix is positive definite.
}
\value{A symmetric matrix of starting values
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
startValues(list(C1,C2))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
