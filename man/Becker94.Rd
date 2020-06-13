\name{Becker94}
\alias{Becker94}
\docType{data}
\title{Five Studies of Ten Correlation Matrices reported by Becker and
  Schram (1994)
}
\description{
  This data set includes five studies of ten correlation matrices
  reported by Becker and Schram (1994). 
}
\usage{data(Becker94)}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 10 correlation matrices. The
	variables are \emph{Math} (math aptitude), \emph{Spatial} (spatial
	ability), and \emph{Verbal} (verbal ability)}
  \item{n}{A vector of sample sizes}
  \item{gender}{\emph{Females} or \emph{Males} samples}
}
}

\source{
  Becker, B. J., & Schram, C. M. (1994). Examining explanatory models through research synthesis. In H. Cooper & L. V. Hedges (Eds.), \emph{The handbook of research synthesis} (pp. 357-381). New York: Russell Sage Foundation.
}

\examples{
\dontrun{
data(Becker94)

#### Fixed-effects model
## First stage analysis
fixed1 <- tssem1(Becker94$data, Becker94$n, method="FEM")
summary(fixed1)

## Prepare a regression model using create.mxMatrix()
## A1 <- create.mxMatrix(c(0,0,0,"0.2*Spatial2Math",
##                         0,0,"0.2*Verbal2Math",0,0), type="Full",
##                       ncol=3, nrow=3, name="A1")
## S1 <- create.mxMatrix(c("0.2*ErrorVarMath",0,0,1,
##                         "0.2*CorBetweenSpatialVerbal",1),
##                       type="Symm", name="S1")

## An alternative method to create a regression model with the lavaan syntax
model <- "## Regression model
          Math ~ Spatial2Math*Spatial + Verbal2Math*Verbal
          ## Error variance of Math
          Math ~~ ErrorVarMath*Math
          ## Variances of Spatial and Verbal fixed at 1.0
          Spatial ~~ 1*Spatial
          Verbal ~~ 1*Verbal
          ## Correlation between Spatial and Verbal
          Spatial ~~ CorBetweenSpatialVerbal*Verbal"

## Display the model
plot(model)

RAM <- lavaan2RAM(model, obs.variables=c("Math", "Spatial", "Verbal"))
RAM

## Second stage analysis
## A1 <- RAM$A
## S1 <- RAM$S
## fixed2 <- tssem2(fixed1, Amatrix=A1, Smatrix=S1, intervals.type="LB")

fixed2 <- tssem2(fixed1, RAM=RAM, intervals.type="LB")
summary(fixed2)

## Display the model with the parameter estimates
plot(fixed2)
    
#### Fixed-effects model: with gender as cluster
## First stage analysis
cluster1 <- tssem1(Becker94$data, Becker94$n, method="FEM", cluster=Becker94$gender)
summary(cluster1)

## Second stage analysis  
cluster2 <- tssem2(cluster1, RAM=RAM, intervals.type="LB")
summary(cluster2)

#### Conventional fixed-effects GLS approach
## First stage analysis
## No random effects
## Replicate Becker's (1992) analysis using 4 studies only
gls1 <- tssem1(Becker92$data[1:4], Becker92$n[1:4], method="REM", RE.type="Zero",
               model.name="Fixed effects GLS Stage 1")
summary(gls1)

## Fixed-effects GLS model: Second stage analysis
gls2 <- tssem2(gls1, RAM=RAM, intervals.type="LB",
               model.name="Fixed effects GLS Stage 2")
summary(gls2)
}
}
\keyword{datasets}
