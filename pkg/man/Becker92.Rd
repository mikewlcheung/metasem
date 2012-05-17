\name{Becker92}
\alias{Becker92}
\docType{data}
\title{Six Studies of Correlation Matrices reported by Becker (1992; 1995)}
\description{
The data set includes six studies of correlation matrices reported by
Becker (1992; 1995).
}
\usage{data(Becker92)}
\format{
}
\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 6 studies of correlation matrices. The
	variables are \emph{Math aptitude}, \emph{Spatial ability} and \emph{Verbal ability}}
  \item{n}{A vector of sample sizes}
}
}
\source{
  Becker, B. J. (1992). Using results from replicated studies to
  estimate linear models. \emph{Journal of Educational Statistics},
  \bold{17(4)}, 341-362. doi:10.3102/10769986017004341

  Becker, B. J. (1995). Corrections to "Using Results from Replicated Studies to Estimate Linear Models." \emph{Journal of Educational and Behavioral Statistics}, \bold{20(1)}, 100-102. doi:10.2307/1165390
}
\references{
}
\examples{
## Show the data
Becker92

#### Fixed-effects model
## First stage analysis
## Replicate Becker's (1992) analysis by using 4 studies only
fixed1 <- tssem1(Becker92$data[1:4], Becker92$n[1:4], method="FEM")
summary(fixed1)

## Prepare a regression model
A1 <- mxMatrix("Full", ncol=3, nrow=3, value=0, free=c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE),
               label=c(NA,NA,NA,"spatial2aptitude",NA,NA,"verbal2aptitude",NA,NA),
               name="A1")
S1 <- mxMatrix("Symm", ncol=3, nrow=3, value=c(1,0,0,1,0.2,1), free=c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE),
               label=c(NA,NA,NA,NA,"cor",NA), name="S1")
Id <- mxMatrix("Iden", ncol=3, nrow=3, name="Id")
## Model implied correlation matrix
impliedR <- mxAlgebra(solve(Id-A1) \%&\% S1, name="impliedR")

## Fixed-effects model: Second stage analysis
fixed2 <- tssem2(fixed1, impliedS=impliedR, matrices=c(A1, Id, S1))
summary(fixed2)

#### Random-effects model
## First stage analysis
## No random effects for off-diagonal elements
random1 <- tssem1(Becker92$data, Becker92$n, method="REM", RE.diag.only=TRUE)
summary(random1)

## Random-effects model: Second stage analysis
random2 <- tssem2(random1, impliedS=impliedR, matrices=c(A1, Id, S1))
summary(random2)
}
\keyword{datasets}
