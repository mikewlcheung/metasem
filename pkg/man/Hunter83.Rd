\name{Hunter83}
\alias{Hunter83}
\docType{data}
\title{Fourteen Studies of Correlation Matrices reported by Hunter (1983)
}
\description{
This data set includes fourteen studies of Correlation Matrices reported by Hunter (1983)
}
\usage{data(Hunter83)}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 14 studies of correlation matrices. The
	variables are \emph{Ability}, \emph{Job knowledge}, \emph{Work
	  sample} and \emph{Supervisor rating}}
  \item{n}{A vector of sample sizes}
}
}
\source{
  Hunter, J. E. (1983). A causal analysis of cognitive ability, job knowledge, job performance, and supervisor ratings. In F. Landy, S. Zedeck, & J. Cleveland (Eds.), \emph{Performance Measurement and Theory} (pp. 257-266). Hillsdale, NJ: Erlbaum.
}

\examples{
data(Hunter83)

#### Fixed-effects model
## First stage analysis
fixed1 <- tssem1(Hunter83$data, Hunter83$n, method="FEM",
                 model.name="TSSEM1 fixed effects model")
summary(fixed1)

#### Second stage analysis
## Model without direct effect from Ability to Supervisor
A1 <- create.mxMatrix(c(0,"0.1*A2J","0.1*A2W",0,0,0,"0.1*J2W","0.1*J2S",0,0,0,"0.1*W2S",0,0,0,0),
                        type="Full", ncol=4, nrow=4)
S1 <- create.mxMatrix(c(1,"0.1*Var_e_J", "0.1*Var_e_W", "0.1*Var_e_S"), type="Diag")

## diag.constraints=TRUE is required as there are mediators  
fixed2 <- tssem2(fixed1, Amatrix=A1, Smatrix=S1, intervals.type="LB", diag.constraints=TRUE,
                 model.name="TSSEM2 fixed effects model")
summary(fixed2)

## Coefficients
coef(fixed2)

## VCOV based on parametric bootstrap
vcov(fixed2)

#### Random-effects model with diagonal elements only
## First stage analysis
random1 <- tssem1(Hunter83$data, Hunter83$n, method="REM", RE.type="Diag", 
                  model.name="TSSEM1 random effects model")
## rerun to remove error code
random1 <- rerun(random1)
summary(random1)

## Second stage analysis
## Model without direct effect from Ability to Supervisor

## diag.constraints=TRUE is required as there are mediators 
random2 <- tssem2(random1, Amatrix=A1, Smatrix=S1, intervals.type="LB", diag.constraints=TRUE,
                  mx.algebras=list( ind=mxAlgebra(A2J*J2S+A2J*J2W*W2S+A2W*W2S, name="ind") ),
                  model.name="TSSEM2 random effects model")
summary(random2)
}
\keyword{datasets}
