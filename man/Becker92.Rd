\name{Becker92}
\alias{Becker92}
\docType{data}
\title{Six Studies of Correlation Matrices reported by Becker (1992; 1995)}
\description{
This data set includes six studies of correlation matrices reported by Becker (1992; 1995).
}
\usage{data(Becker92)}

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

\examples{
\dontrun{
data(Becker92)

#### Fixed-effects model
## First stage analysis
## Replicate Becker's (1992) analysis using 4 studies only
fixed1 <- tssem1(Becker92$data[1:4], Becker92$n[1:4], method="FEM")
summary(fixed1)

## Prepare a regression model using create.mxMatrix()
A1 <- create.mxMatrix(c(0,0,0,"0.2*Spatial2Math",
                        0,0,"0.2*Verbal2Math",0,0), type="Full",
                        ncol=3, nrow=3, as.mxMatrix=FALSE)

## This step is not necessary but it is useful for inspecting the model.
dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- c("Math","Spatial","Verbal") 

## Display A1
A1

S1 <- create.mxMatrix(c("0.2*ErrorVarMath",0,0,1,"0.2*CorSpatialVerbal",1),
                        type="Symm", as.mxMatrix=FALSE)

## This step is not necessary but it is useful for inspecting the model.
dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- c("Math","Spatial","Verbal") 

## Display S1
S1

################################################################################
## Alternative model specification in lavaan model syntax
model <- "## Regression paths
          Math~Spatial2Math*Spatial
          Math~Verbal2Math*Verbal
          Spatial~~CorSpatialVerbal*Verbal
          ## Fix the variances of Spatial and Verbal at 1
          Spatial~~1*Spatial
          Verbal~~1*Verbal
          ## Label the error variance of Math
          Math~~ErrorVarMath*Math+start(0.2)*Math"

## Display the model
plot(model)

RAM <- lavaan2RAM(model, obs.variables=c("Math","Spatial","Verbal"))
RAM

A1 <- RAM$A
S1 <- RAM$S
################################################################################

## Fixed-effects model: Second stage analysis
## Two equivalent versions to calculate the R2 and its 95\% LBCI
fixed2 <- tssem2(fixed1, Amatrix=A1, Smatrix=S1, intervals.type="LB",
       mx.algebras=list(R1=mxAlgebra(Spatial2Math^2+Verbal2Math^2
                        +2*CorSpatialVerbal*Spatial2Math*Verbal2Math, name="R1"),
                        R2=mxAlgebra(One-Smatrix[1,1], name="R2"),
                        One=mxMatrix("Iden", ncol=1, nrow=1, name="One")))
summary(fixed2)

## Display the model with the parameter estimates
plot(fixed2)

#### Random-effects model
## First stage analysis
## No random effects for off-diagonal elements
random1 <- tssem1(Becker92$data, Becker92$n, method="REM", RE.type="Diag")
summary(random1)

## Random-effects model: Second stage analysis
random2 <- tssem2(random1, Amatrix=A1, Smatrix=S1, intervals.type="z")
summary(random2)

## Display the model with the parameter estimates
plot(random2, color="yellow")

#### Conventional fixed-effects GLS approach
## First stage analysis
## No random effects
## Replicate Becker's (1992) analysis using 4 studies only
gls1 <- tssem1(Becker92$data[1:4], Becker92$n[1:4], method="REM", RE.type="Zero",
               model.name="Fixed effects GLS Stage 1")
summary(gls1)

## Fixed-effects GLS model: Second stage analysis
gls2 <- tssem2(gls1, Amatrix=A1, Smatrix=S1, intervals.type="z",
               model.name="Fixed effects GLS Stage 2")
summary(gls2)
}
}
\keyword{datasets}
