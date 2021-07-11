\name{Hunter83}
\alias{Hunter83}
\docType{data}
\title{Fourteen Studies of Correlation Matrices reported by Hunter (1983)
}
\description{
This dataset includes fourteen studies of Correlation Matrices reported by Hunter (1983)
}
\usage{data(Hunter83)}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 14 studies of correlation matrices. The
	variables are \emph{Ability}, \emph{Job knowledge},
	\emph{Work sample} and \emph{Supervisor rating}}
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
## A1 <- create.mxMatrix(c(0,"0.1*A2J","0.1*A2W",0,0,0,"0.1*J2W","0.1*J2S",
##                         0,0,0,"0.1*W2S",0,0,0,0),
##                         type="Full", ncol=4, nrow=4, as.mxMatrix=FALSE)

## ## This step is not necessary but it is useful for inspecting the model.
## dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- c("Ability","Job","Work","Supervisor") 
## A1

## S1 <- create.mxMatrix(c(1,"0.1*Var_e_J", "0.1*Var_e_W", "0.1*Var_e_S"),
##                       type="Diag", as.mxMatrix=FALSE)
## dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- c("Ability","Job","Work","Supervisor") 
## S1

################################################################################
## Model specification in lavaan model syntax
## The "ind" effect can be defined within the syntax
model1 <- "## Regression paths
           Job_knowledge ~ A2J*Ability
           Work_sample ~ A2W*Ability + J2W*Job_knowledge
           Supervisor ~ J2S*Job_knowledge + W2S*Work_sample

           ## Fix the variance of Ability at 1
           Ability ~~ 1*Ability

           ## Label the error variances of the dependent variables
           Job_knowledge ~~ VarE_J*Job_knowledge
           Work_sample ~~ VarE_W*Work_sample
           Supervisor ~~ VarE_S*Supervisor

           ## Define an indirect effect
           ind := A2J*J2S+A2J*J2W*W2S+A2W*W2S"

## Display the model
plot(model1, layout="spring", sizeMan=10)

RAM1 <- lavaan2RAM(model1, obs.variables=c("Ability","Job_knowledge",
                   "Work_sample","Supervisor"))
RAM1

################################################################################
fixed2 <- tssem2(fixed1, RAM=RAM1, intervals.type="z",
                 diag.constraints=FALSE,
                 model.name="TSSEM2 fixed effects model")
summary(fixed2)

## Display the model with the parameter estimates
plot(fixed2, layout="spring")

## Coefficients
coef(fixed2)

## VCOV based on parametric bootstrap
vcov(fixed2)

#### Random-effects model with diagonal elements only
## First stage analysis
random1 <- tssem1(Hunter83$data, Hunter83$n, method="REM", RE.type="Diag", 
                  acov="individual", model.name="TSSEM1 random effects model")
summary(random1)

model2 <- "## Regression paths
           Job_knowledge ~ A2J*Ability
           Work_sample ~ A2W*Ability + J2W*Job_knowledge
           Supervisor ~ J2S*Job_knowledge + W2S*Work_sample

           ## Fix the variance of Ability at 1
           Ability ~~ 1*Ability

           ## Label the error variances of the dependent variables
           Job_knowledge ~~ VarE_J*Job_knowledge
           Work_sample ~~ VarE_W*Work_sample
           Supervisor ~~ VarE_S*Supervisor"

RAM2 <- lavaan2RAM(model2, obs.variables=c("Ability","Job_knowledge",
                   "Work_sample","Supervisor"))
RAM2

## Second stage analysis
## Model without direct effect from Ability to Supervisor

## The "ind" effect is defined in tssem2().
random2 <- tssem2(random1, RAM=RAM2, intervals.type="LB",
                  diag.constraints=FALSE,
                  mx.algebras=
                      list(ind=mxAlgebra(A2J*J2S+A2J*J2W*W2S+A2W*W2S, name="ind")),
                  model.name="TSSEM2 random effects model")
                  
summary(random2)

## Display the model with the parameter estimates
plot(random2, layout="spring")
}
}
\keyword{datasets}
