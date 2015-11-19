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
\dontrun{
data(Hunter83)

#### Fixed-effects model
## First stage analysis
fixed1 <- tssem1(Hunter83$data, Hunter83$n, method="FEM",
                 model.name="TSSEM1 fixed effects model")
summary(fixed1)

#### Second stage analysis
## Model without direct effect from Ability to Supervisor
A1 <- create.mxMatrix(c(0,"0.1*A2J","0.1*A2W",0,0,0,"0.1*J2W","0.1*J2S",
                        0,0,0,"0.1*W2S",0,0,0,0),
                        type="Full", ncol=4, nrow=4, as.mxMatrix=FALSE)

## This step is not necessary but it is useful for inspecting the model.
dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- c("Ability","Job","Work","Supervisor") 
A1

S1 <- create.mxMatrix(c(1,"0.1*Var_e_J", "0.1*Var_e_W", "0.1*Var_e_S"),
                      type="Diag", as.mxMatrix=FALSE)
dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- c("Ability","Job","Work","Supervisor") 
S1

################################################################################
## Alternative model specification in lavaan model syntax
model <- "## Regression paths
          Job~A2J*Ability
          Work~A2W*Ability + J2W*Job
          Supervisor~J2S*Job + W2S*Work
          ## Fix the variance of Ability
          Ability~~1*Ability
          ## Label the error variances of dependent variables
          Job~~Var_e_J*Job
          Work~~Var_e_W*Work
          Supervisor~~Var_e_S*Supervisor"

RAM <- lavaan2RAM(model, obs.variables=c("Ability","Job","Work","Supervisor"))
RAM

A1 <- RAM$A
S1 <- RAM$S
################################################################################

## diag.constraints=TRUE is required as there are mediators  
fixed2 <- tssem2(fixed1, Amatrix=A1, Smatrix=S1, intervals.type="LB",
                 diag.constraints=FALSE,
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
summary(random1)

## Second stage analysis
## Model without direct effect from Ability to Supervisor

## diag.constraints=TRUE is required as there are mediators 
random2 <- tssem2(random1, Amatrix=A1, Smatrix=S1, intervals.type="LB",
                  diag.constraints=FALSE,
                  mx.algebras=
                  list( ind=mxAlgebra(A2J*J2S+A2J*J2W*W2S+A2W*W2S, name="ind") ),
                  model.name="TSSEM2 random effects model")
summary(random2)

## Load the library
library("semPlot")

## Convert the model to semPlotModel object
my.plot <- meta2semPlot(random2)

## Plot the model with labels
semPaths(my.plot, whatLabels="path", nCharEdges=10, nCharNodes=10, color="red")

## Plot the parameter estimates
semPaths(my.plot, whatLabels="est", nCharNodes=10, color="green")
}
}
\keyword{datasets}
