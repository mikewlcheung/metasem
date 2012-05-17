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
issp89

#### Fixed-effects model: Stage 1 analysis
fixed1 <- tssem1(issp89$data, issp89$n, method="FEM")
summary(fixed1)
  
## Prepare a model implied matrix
## Factor correlation matrix
Phi <- mxMatrix("Stand", ncol=3, nrow=3, value=.2, free=TRUE, name="Phi")

## Factor loadings
Lambda <- as.mxMatrix( matrix(c(".3*a",".3*b",".3*c",rep(0,9),
                                ".3*d",".3*e",".3*f",".3*g",rep(0,9),
                                ".3*h",".3*i"), ncol=3, nrow=9), name="Lambda")

impliedR <- mxAlgebra(Lambda \%&\% Phi, name="impliedR")
  
#### Fixed-effects model: Stage 2 analysis
fixed2 <- tssem2(fixed1, impliedS=impliedR, matrices=c(Phi, Lambda))
summary(fixed2)
  
#### Random-effects model: Stage 1 analysis
## Non-convergent solution
% random1 <- tssem1(issp89$data, issp89$n, method="REM", RE.diag.only=TRUE)
% summary(random1)
  
% #### Random-effects model: Stage 2 analysis
% random2 <- tssem2(random1, impliedS=impliedR, matrices=c(Phi, Lambda))
% summary(random2) 
}
\keyword{datasets}
