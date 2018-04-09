\name{Cooke16}
\alias{Cooke16}
\docType{data}
\title{Correlation Matrices from Cooke et al. (2016)}
\description{
  The data set includes correlation matrices on using the theory of planned
  behavior to predict alcohol consumption reported by Cooke et al. (2016). 
}
\usage{data(Cooke16)}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of correlation matrices. The variables are
	\emph{SN} (subjective norm), \emph{ATT} (attitude), \emph{PBC}
	(perceived behavior control), \emph{BI} (behavioral intention), and
	\emph{BEH} (behavior)}
  \item{n}{A vector of sample sizes}
}
}
\source{
  Cooke, R., Dahdah, M., Norman, P., & French, D. P. (2016). How well does the theory of planned behaviour predict alcohol consumption? A systematic review and meta-analysis. \emph{Health Psychology Review}, \bold{10}(2), 148-167.
}
\references{
  Cheung, M. W.-L., & Hong, R. Y. (2017). Applications of meta-analytic structural equation modeling in health psychology: Examples, issues, and recommendations. \emph{Health Psychology Review}, \bold{11}, 265-279. 
}

\examples{
\dontrun{
## Check whether the correlation matrices are valid (positive definite)
Cooke16$data[is.pd(Cooke16$data)==FALSE]

## Since the correlation matrix in Study 3 is not positive definite,
## we exclude it in the following analyses
my.data <- Cooke16$data[-3]
my.n <- Cooke16$n[-3]

## Show the no. of studies per correlation
pattern.na(my.data, show.na = FALSE)

## Show the total sample sizes per correlation
pattern.n(my.data, my.n)

## Stage 1 analysis
## Random-effects model
random1 <- tssem1(my.data, my.n, method="REM", RE.type="Diag")
summary(random1)

A1 <- create.mxMatrix(c(0,0,0,0,0,
                        0,0,0,0,0,
                        0,0,0,0,0,
                        "0.2*SN2BI","0.2*ATT2BI","0.2*PBC2BI",0,0,
                        0,0,"0.2*PBC2BEH","0.2*BI2BEH",0),
                        type="Full", ncol=5, nrow=5,
                        byrow=TRUE, as.mxMatrix=FALSE)

## This step is not necessary but it is useful for inspecting the model.
dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- colnames(Cooke16$data[[1]])

## Display A1
A1

S1 <- create.mxMatrix(c(1,
                        "0.1*ATT_SN", 1,
                        "0.1*PBC_SN", "0.1*PBC_ATT", 1,
                        0, 0, 0, "0.5*VarBI",
                        0, 0, 0, 0, "0.5*VarBEH"),
                      type = "Symm", ncol=5, nrow=5,
                      byrow=TRUE, as.mxMatrix=FALSE)

dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- colnames(Cooke16$data[[1]])
S1

## Stage 2 analysis
random2 <- tssem2(random1, Amatrix=A1, Smatrix=S1, diag.constraints=FALSE,
                  intervals.type="LB")
summary(random2)

## Display the model
plot(random2, what="path")    
    
## Display the model with the parameter estimates
plot(random2, color="yellow")
}
}
\keyword{datasets}
