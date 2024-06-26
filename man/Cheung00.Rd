\name{Cheung00}
\alias{Cheung00}
\docType{data}
\title{Fifty Studies of Correlation Matrices used in Cheung and Chan (2000)}
\description{
This data set includes fifty studies of correlation matrices on the theory
of planned theory reported by Cheung and Chan (2000).
}
\usage{data(Cheung00)}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 50 studies of correlation matrices. The
	variables are the attitude toward behavior \emph{att}, subjective norm \emph{sn},
	behavioral intention \emph{bi}, and behavior \emph{beh}}
  \item{n}{A vector of sample sizes}
}
}
\source{
  Cheung, S.-F., & Chan, D. K.-S. (2000). The role of perceived behavioral control in predicting human behavior: A meta-analytic review of studies on the theory of planned behavior. \emph{Unpublished manuscript}, Chinese University of Hong Kong.
}
\references{
  Cheung, M.W.-L., & Cheung, S.-F. (2016). Random-effects models for
  meta-analytic structural equation modeling: Review, issues, and
  illustrations. \emph{Research Synthesis Methods}, \bold{7}, 140-155.
}
\note{These studies were extracted from the original data set for
  illustration purpose. Some samples contained two or more correlation
  matrices, and only one of them was arbitrarily selected to avoid the
  problem of dependence. Moreover, studies with less than 3 correlation
  coefficients were also excluded.
  }
\examples{
\donttest{
data(Cheung00)

## Variable labels
labels <- colnames(Cheung00$data[[1]])

## Full mediation model
S <- create.mxMatrix(c("1", 
                       ".2*cov_att_sn", "1", 
                       0, 0, ".2*e_bi", 
                       0, 0, 0, ".2*e_beh"), 
                     type="Symm", as.mxMatrix=FALSE, byrow=TRUE)
dimnames(S) <- list(labels, labels)
S

A <- matrix(c("0","0","0","0",
              "0","0","0","0",
              ".2*att2bi", ".2*sn2bi", "0", "0",
              "0", "0", ".2*bi2beh", "0"),
            byrow=TRUE, 4, 4)
dimnames(A) <- list(labels, labels)
A

#### Random-effects model

## Stage 1 analysis
random_1 <- tssem1(Cheung00$data, Cheung00$n, method="REM", RE.type="Diag",
                   acov="weighted")
summary(random_1)

## Stage 2 analysis
random_2 <- tssem2(random_1, Amatrix=A, Smatrix=S, intervals.type="LB",
                   diag.constraints=TRUE)
summary(random_2)

## Display the model
plot(random_2, what="path")

## Display the model with the parameter estimates
plot(random_2, color="yellow")

## Load the library
library("semPlot")
}
}
\keyword{datasets}
