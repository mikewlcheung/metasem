\name{issp05}
\alias{issp05}
\docType{data}
\title{Data Set from ISSP (2005) 
}
\description{Thirty-two covariance matrices on work-related attitudes
  were extracted from the International Social Survey Programme 2005:
  Work Orientation III (ISSP, 2005). Seven variables were selected for
  demonstration purposes. They were grouped into three constructs:
  \emph{Importance of Job Prospects} measured by job security (JP1),
  high income (JP2), and opportunity for advancement (JP3); \emph{Importance
  of Job Autonomy} measured by work independently (JA1) and decide time
  of work (JA2); and \emph{Importance of Contributions to Society} measured by help other
  people (CS1) and a job useful to society (CS2).
}
\usage{data(issp05)}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 32 covariance matrices }
  \item{n}{A vector of sample sizes}
}
}
\source{ISSP Research Group (2007): \emph{International Social Survey Programme 2005: Work Orientation III (ISSP 2005)}. GESIS Data Archive, Cologne. ZA4350 Data file Version 1.0.0, doi:10.4232/1.4350}


\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{issp89}}
}

\examples{

\dontrun{
data(issp05)

#### Fixed-effects TSSEM
fixed1 <- tssem1(issp05$data, issp05$n, method="FEM")
summary(fixed1)

## Prepare a model implied matrix
## Factor correlation matrix
Phi <- create.mxMatrix( c("0.3*corf2f1","0.3*corf3f1","0.3*corf3f2"),
                        type="Stand", as.mxMatrix=FALSE )

## Error variances
Psi <- create.mxMatrix( paste("0.2*e", 1:7, sep=""), type="Diag",
                        as.mxMatrix=FALSE )

## Create Smatrix
S1 <- bdiagMat(list(Psi, Phi))
## dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- c(paste("x",1:7,sep=""),
##                                             paste("f",1:3,sep=""))
## S1
S1 <- as.mxMatrix(S1)

## Factor loadings
Lambda <- create.mxMatrix( c(".3*f1x1",".3*f1x2",".3*f1x3",rep(0,7),
                             ".3*f2x4",".3*f2x5",rep(0,7),".3*f3x6",".3*f3x7"),
                           type="Full", ncol=3, nrow=7, as.mxMatrix=FALSE )
Zero1 <- matrix(0, nrow=7, ncol=7)
Zero2 <- matrix(0, nrow=3, ncol=10)

## Create Amatrix
A1 <- rbind( cbind(Zero1, Lambda),
             Zero2 )
## dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- c(paste("x",1:7,sep=""),
##                                             paste("f",1:3,sep=""))
## A1
A1 <- as.mxMatrix(A1)

## Create Fmatrix
F1 <- create.Fmatrix(c(rep(1,7), rep(0,3)))

#### Fixed-effects model: Stage 2 analysis
fixed2 <- tssem2(fixed1, Amatrix=A1, Smatrix=S1, Fmatrix=F1,
                 intervals.type="LB")
summary(fixed2)
}
}
\keyword{datasets}
