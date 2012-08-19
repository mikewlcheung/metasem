\name{Cheung09}
\alias{Cheung09}
\docType{data}
\title{A Data Set from TSSEM User's Guide Version 1.11 by Cheung (2009) 
}
\description{Four studies were selected from the data set used by Cheung
  and Chan (2005; 2009). Some variables were randomly deleted to
  illustrate the analysis with missing data.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(Cheung09)}
\format{
}
\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 4 studies of correlation matrices }
  \item{n}{A vector of sample sizes}
}
}
\source{Cheung, M. W. L. (2009). TSSEM: A LISREL syntax generator for two-stage structural equation modeling (Version 1.11) [Computer software]. Retrieved from \url{http://courses.nus.edu.sg/course/psycwlm/internet/tssem.zip}.
}
\references{
  Cheung, M. W. L., & Chan, W. (2005). Meta-analytic structural equation modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10}, 40-64.

Cheung, M. W. L., & Chan, W. (2009). A two-stage approach to synthesizing covariance matrices in meta-analytic structural equation modeling. \emph{Structural Equation Modeling}, \bold{16}, 28-53.
}
\examples{
data(Cheung09)

#### Fixed-effects model: Stage 1 analysis
fixed1 <- tssem1(Cheung09$data, Cheung09$n, method="FEM")
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
