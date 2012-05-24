\name{issp89}
\alias{issp89}
\docType{data}
\title{Data Set from Cheung and Chan (2005; 2009) 
}
\description{Eleven covariance matrices on work-related attitudes were
  extracted from the Inter-University Consortium for Political and
  Social Research (1989). Nine variables were selected by Cheung and
  Chan (2005; 2009) for demonstration purposes. They were grouped into
  three constructs: \emph{Job Prospects} measured by job security
  (JP1), income (JP2), and advancement opportunity (JP3); \emph{Job Nature}
  measured by interesting job (JN1), independent work (JN2), help
  other people (JN3), and useful to society (JN4); and \emph{Time Demand} measured by flexible working hours (TD1) and lots of leisure time (TD2).
}
\usage{data(issp89)}
\format{
}
\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 11 studies of covariance matrices }
  \item{n}{A vector of sample sizes}
}
}
\source{Inter-University Consortium for Political and Social Research. (1989). \emph{International Social Survey Program: Work orientation}. Ann Arbor, MI: Author.}

\references{
  Cheung, M. W. L., & Chan, W. (2005). Meta-analytic structural equation modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10}, 40-64.

Cheung, M. W. L., & Chan, W. (2009). A two-stage approach to synthesizing covariance matrices in meta-analytic structural equation modeling. \emph{Structural Equation Modeling}, \bold{16}, 28-53.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{issp05}}
}
\examples{
data(issp89)

#### Fixed-effects model: Stage 1 analysis
fixed1 <- tssem1(issp89$data, issp89$n, method="FEM")
summary(fixed1)
  
## Prepare a model implied matrix
## Factor correlation matrix
Phi <- create.mxMatrix( c("0.3*corf2f1","0.3*corf3f1","0.3*corf3f2"),
                        type="Stand", as.mxMatrix=FALSE )
## Error variances
Psi <- create.mxMatrix( paste("0.2*e", 1:9, sep=""), type="Diag",
                        as.mxMatrix=FALSE )

## Create Smatrix
S1 <- bdiagMat(list(Psi, Phi))
## dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- c(paste("x",1:9,sep=""), paste("f",1:3,sep=""))
## S1
S1 <- as.mxMatrix(S1)

## Factor loadings
Lambda <- create.mxMatrix( c(".3*f1x1",".3*f1x2",".3*f1x3",rep(0,9),
                             ".3*f2x4",".3*f2x5",".3*f2x6",".3*f2x7",
                             rep(0,9),".3*f3x8",".3*f3x9"), type="Full",
                             ncol=3, nrow=9, as.mxMatrix=FALSE )
Zero1 <- matrix(0, nrow=9, ncol=9)
Zero2 <- matrix(0, nrow=3, ncol=12)

## Create Amatrix
A1 <- rbind( cbind(Zero1, Lambda),
             Zero2 )
## dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- c(paste("x",1:9,sep=""), paste("f",1:3,sep=""))
## A1
A1 <- as.mxMatrix(A1)

## Create Fmatrix
F1 <- create.Fmatrix(c(rep(1,9), rep(0,3)))
  
#### Fixed-effects model: Stage 2 analysis
fixed2 <- tssem2(fixed1, Amatrix=A1, Smatrix=S1, Fmatrix=F1,
                 diag.constraints=TRUE, intervals.type="LB")
summary(fixed2)
}
\keyword{datasets}
