\name{impliedSigma}
\alias{impliedSigma}
\title{Create the Model Implied Correlation or Covariance Matrix
}
\description{It creates the model implied correlation or covariance
  matrix based on the RAM model specification.
}
\usage{
impliedSigma(Amatrix=NULL, Smatrix=NULL, Fmatrix=NULL, labels=NULL,
             cor.analysis=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Amatrix}{An asymmetric matrix in the RAM specification with
	\code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be
	converted into \code{\link[OpenMx]{MxMatrix-class}} by the
	\code{as.mxMatrix} function.
}
  \item{Smatrix}{A symmetric matrix in the RAM specification with
	\code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be
	converted into \code{\link[OpenMx]{MxMatrix-class}} by the
	\code{as.mxMatrix} function.
}
  \item{Fmatrix}{A filter matrix in the RAM specification with
	\code{\link[OpenMx]{MxMatrix-class}}. If it is \code{NULL} (the default), an
	identity matrix with the same dimensions of \code{Smatrix} will be
	created. If it is a matrix, it will be converted into
	\code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix}
	function. It is not required when there is no latent variable.
  }
  \item{labels}{A character vectors of labels of observed and latent
	variables with the same dimensions as that in the \code{Amatrix} and
	\code{Smatrix}. The labels of the observed variables are
	automatically extracted by the \code{Fmatrix}.
  }
  \item{cor.analysis}{Logical. The output is either the model implied
	correlation matrix or the covariance matrix.
  }
  \item{\dots}{Not used.}
}
\details{This function can be used to generate the model implied
  correlation matrix for the standardized parameters with the
  \code{cor.analysis=TRUE} argument. Suppose we want to calculate the
  population correlation matrix for a mediation model with x, m, and
  y. We only need to specify the population path coefficients among x,
  m, and y in the \code{Amatrix}. We do not need to specify the
  population error variances of m and y. We treat the error variances as
  unknown parameters by giving them starting values in the
  \code{Smatrix} matrix. When the covariance matrix is requested by
  specifying \code{cor.analysis=FALSE}, it simply calculates the
  population model covariance matrix by treating the values in
  \code{Smatrix} as the population values.
}
\note{It is important to ensure that all the population values in
  \code{Amatrix} must be set as fixed parameters; otherwise, these
  values may be altered with the \code{cor.analysis=TRUE} argument. When
  there is an error or warning message about the status code, there is a
  high chance that some of the values in \code{Amatrix} are incorrectly
  set as free parameters.
}
  
\value{A list of RAM matrices, the model implied correlation or
  covariance matrix of the observed variables (\code{SigmaObs}), of both
  observed and latent variables (\code{SigmaAll}), the minimum fit
  (\code{minFit}) which should be zero, and the status code of the
  optimization (\code{status}) which should also be zero. The last
  object is \code{mx.fit} which is the output after running the
  model. It can be used for diagnosis.
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\examples{
## A simple mediation model
## All are population parameters in the A matrix
A1 <- matrix(c(0, 0, 0,
               0.3, 0, 0,
               0.2, 0.4, 0), nrow=3, ncol=3, byrow=TRUE,
             dimnames=list(c("x", "m", "y"), c("x", "m", "y")))
A1             

## Variance of x is fixed at 1 while the other variances are free.
S1 <- matrix(c(1, 0, 0,
               0, "0.1*ErrVarM",0,
               0, 0, "0.1*ErrVarY"), nrow=3, ncol=3,
             dimnames=list(c("x", "m", "y"), c("x", "m", "y")))
S1

impliedSigma(Amatrix=A1, Smatrix=S1, labels=c("x","m","y"))

## A CFA model
A2 <- matrix(c(0, 0, 0, 0.3,
               0, 0, 0, 0.4,
               0, 0, 0, 0.5,
               0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE,
            dimnames=list(c("x1", "x2", "x3", "f"), c("x1", "x2", "x3", "f")))
A2

## Variance of f is fixed at 1 while the other variances are free.
S2 <- matrix(c("0.1*Err1", 0, 0, 0,
                0, "0.1*Err2", 0, 0,
                0, 0, "0.1*Err3", 0,
                0, 0, 0, 1), nrow=4, ncol=4,
            dimnames=list(c("x1", "x2", "x3", "f"), c("x1", "x2", "x3", "f")))
S2

F2 <- create.Fmatrix(c(1,1,1,0), as.mxMatrix=FALSE)
F2

impliedSigma(Amatrix=A2, Smatrix=S2, Fmatrix=F2, labels=c("x1", "x2", "x3", "f"))
}
\keyword{utilities}
