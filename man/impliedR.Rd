\name{impliedR}
\alias{impliedR}
\alias{rimpliedR}
\title{Create or Generate the Model Implied Correlation or Covariance Matrices
}
\description{It creates or generates the model implied correlation or covariance
  matrices based on the RAM model specification.
}
\usage{
impliedR(Amatrix, Smatrix, Fmatrix, corr=TRUE, labels, ...)
rimpliedR(Amatrix, Smatrix, Fmatrix, AmatrixSD, k=1, corr=TRUE,
          nonPD.pop=c("replace", "nearPD", "accept"))
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
	\code{\link[OpenMx]{MxMatrix-class}}. If it is missing, an
	identity matrix with the same dimensions of \code{Smatrix} will be
	created, which means that all variables are observed. If it is a matrix, it will be converted into
	\code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix}
	function. It is not required when there is no latent variable.
  }
  \item{AmatrixSD}{Standard deviations (SD) of the elements in the
	\code{Amatrix}. If it is missing, a matrix of zero is created.}
  \item{k}{Number of studies.}  
  \item{corr}{Logical. The output is either the model implied
	correlation matrix or the covariance matrix.
  }
  \item{labels}{A character vector of the observed and latent
	variables with the same dimensions as that in the \code{Amatrix} and
	\code{Smatrix}. 
  }
  \item{nonPD.pop}{If it is \code{replace}, generated non-positive
	definite matrices are replaced by generated new ones which are
	positive definite. If it is \code{nearPD}, they are replaced by
	nearly positive definite matrices by calling
	\code{Matrix::nearPD()}. If it is \code{accept}, they are accepted.}
  \item{\dots}{Not used.}
}
\details{This function can be used to generate the model implied
  correlation matrix for the standardized parameters with the
  \code{corr=TRUE} argument. Suppose we want to calculate the
  population correlation matrix for a mediation model with x, m, and
  y. We only need to specify the population path coefficients among x,
  m, and y in the \code{Amatrix}. We do not need to specify the
  population error variances of m and y. We treat the error variances as
  unknown parameters by giving them starting values in the
  \code{Smatrix} matrix. When the covariance matrix is requested by
  specifying \code{corr=FALSE}, it simply calculates the
  population model covariance matrix by treating the values in
  \code{Smatrix} as the population values.
}
\note{It is important to ensure that all the population values in
  \code{Amatrix} must be set as fixed parameters; otherwise, these
  values may be altered with the \code{corr=TRUE} argument. When
  there is an error or warning message about the status code, there is a
  high chance that some of the values in \code{Amatrix} are incorrectly
  set as free parameters.
}
  
\value{A list of RAM matrices, the model implied correlation or
  covariance matrix of the observed variables (\code{SigmaObs}), of both
  observed and latent variables (\code{SigmaAll}), the minimum fit
  (\code{minFit}) which should be zero, and the status code of the
  optimization (\code{status}) which should also be zero when the
  optimization is fine. The last object is \code{mx.fit} which is the output after running the
  model. It can be used in the diagnosis.
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\examples{
set.seed(100)

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

impliedR(Amatrix=A1, Smatrix=S1)

## SD of A1
A1SD <- matrix(c(0, 0, 0,
                 0.1, 0, 0,
                 0.1, 0.1, 0), nrow=3, ncol=3, byrow=TRUE,
               dimnames=list(c("x", "m", "y"), c("x", "m", "y")))
A1SD

rimpliedR(Amatrix=A1, Smatrix=S1, AmatrixSD=A1SD, k=2)

## A CFA model
A2 <- matrix(c(0, 0, 0, 0.3,
               0, 0, 0, 0.4,
               0, 0, 0, 0.5,
               0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE,
             dimnames=list(c("x1", "x2", "x3", "f"),
                           c("x1", "x2", "x3", "f")))
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

impliedR(Amatrix=A2, Smatrix=S2, Fmatrix=F2)

## SD of A2
A2SD <- matrix(c(0, 0, 0, 0.1,
                 0, 0, 0, 0.1,
                 0, 0, 0, 0.1,
                 0, 0, 0, 0), nrow=4, ncol=4, byrow=TRUE,
               dimnames=list(c("x1", "x2", "x3", "f"),
                             c("x1", "x2", "x3", "f")))               
A2SD

rimpliedR(Amatrix=A2, Smatrix=S2, Fmatrix=F2, AmatrixSD=A2SD, k=2)
}
\keyword{utilities}
