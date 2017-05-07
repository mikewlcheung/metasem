\name{rCor}
\alias{rCor}
\alias{rCorPop}
\alias{rCorSam}
\title{Generate Sample/Population Correlation/Covariance Matrices
}
\description{It generates random sample or population correlation or
  covariance matrices. \code{rCor()} is a wrapper to call
  \code{rCorPop()} and then \code{rCorSam()}.
}
\usage{
rCor(Sigma, V, n, corr=TRUE, raw.data=FALSE,
     nonPD.pop=c("replace", "nearPD", "accept"),
     nonPD.sam=c("stop", "nearPD"))
rCorPop(Sigma, V, k, corr=TRUE, 
        nonPD.pop=c("replace", "nearPD", "accept"))
rCorSam(Sigma, n, corr=TRUE, raw.data=FALSE, 
        nonPD.sam=c("stop", "nearPD"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sigma}{A list of population correlation/covariance matrices or a
  single matrix}
  \item{V}{A variance-covariance matrix of Sigma.}
  \item{n}{A vector or a single sample sizes.}
  \item{corr}{Logical. Whether to generate correlation or covariance matrices.}
  \item{raw.data}{Logical. Whether correlation/covariance matrices are
	generated via raw.data or directly from a Wishart distribution.}
  \item{nonPD.pop}{If it is \code{replace}, generated non-positive
	definite matrices are replaced by generated new ones which are
	positive definite. If it is \code{nearPD}, they are replaced by
	nearly positive definite matrices by calling
	\code{Matrix::nearPD()}. If it is \code{accept}, they are accepted.}
  \item{nonPD.sam}{If it is \code{stop}, the program stops when the
	inputs in the \code{rCorSam} are non-positive definite. If it is
	\code{nearPD}, they are replaced by nearly positive definite matrices by calling
	\code{Matrix::nearPD()}.}
  \item{k}{A vector or a single number of studies.}
}

\value{An object of the generated population/sample
  correlation/covariance matrices.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\examples{
Sigma <- matrix(c(1, .2, .3,
                  .2, 1, .4,
                  .3, .4, 1), ncol=3, nrow=3)
V <- diag(c(.1, .1, .1))

## Generate two population correlation matrices
Pop.corr <- rCorPop(Sigma, V, k=2)
Pop.corr
                  
summary(Pop.corr)
                  
## Generate two sample correlation matrices
rCorSam(Sigma=Pop.corr, n=c(10, 10))
                 
## The above code is the same as the following one
rCor(Sigma, V, n=c(10, 10))                  
}
\keyword{ utilities }
