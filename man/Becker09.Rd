\name{Becker09}
\alias{Becker09}
\docType{data}
\title{Ten Studies of Correlation Matrices used by Becker (2009)
}
\description{This dataset includes ten studies on the relationships between CSAI
  subscales and sports behavior. The original data were used in Craft et
  al. (2003), whereas a subset of them was illustrated in Becker (2009).
}
\usage{data("Becker09")}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 4x4 correlation matrices. The
	variables are \emph{Performance}, \emph{Cognitive}, \emph{Somatic},
	and \emph{Self_confidence}}
  \item{n}{A vector of sample sizes}
  \item{Type_of_sport}{Samples based on \emph{Individual} or \emph{Team}}
}
}

\source{
Craft, L. L., Magyar, T. M., Becker, B. J., & Feltz, D. L. (2003). The relationship between the Competitive State Anxiety Inventory-2 and sport performance: a meta-analysis. \emph{Journal of Sport and Exercise Psychology}, \bold{25(1)}, 44-65.
}
\references{
Becker, B. J. (2009). Model-based meta-analysis. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 377-395). New York: Russell Sage Foundation.
}
\examples{
\dontrun{
data(Becker09)

#### Fixed-effects model
## First stage analysis
fixed1 <- tssem1(Becker09$data, Becker09$n, method="FEM")
summary(fixed1)

## Prepare a regression model using create.mxMatrix()
A1 <- create.mxMatrix(c(0, "0.1*Cog2Per", "0.1*SO2Per", "0.1*SC2Per",
                        0, 0, 0, 0,
                        0, 0, 0, 0,
                        0, "0.1*Cog2SC", "0.1*SO2SC",0),
                      type="Full", byrow=TRUE, ncol=4, nrow=4,
                      as.mxMatrix=FALSE)

## This step is not necessary but it is useful for inspecting the model.
dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- c("Performance", "Cognitive",
                                            "Somatic", "Self_confidence")

## Display A1
A1

S1 <- create.mxMatrix(c("0.1*var_Per",
                        0, 1,
                        0, "0.1*cor", 1,
                        0, 0, 0, "0.1*var_SC"), byrow=TRUE, type="Symm",
                      as.mxMatrix=FALSE)

## This step is not necessary but it is useful for inspecting the model.
dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- c("Performance", "Cognitive",
                                            "Somatic", "Self_confidence")

## Display S1
S1

################################################################################
## Alternative model specification in lavaan model syntax
model <- "## Regression paths
          Performance ~ Cog2Per*Cognitive + SO2Per*Somatic + SC2Per*Self_confidence
          Self_confidence ~ Cog2SC*Cognitive + SO2SC*Somatic
          ## Fix the variances of Cog and SO at 1
          Cognitive ~~ 1*Cognitive
          Somatic ~~ 1*Somatic
          ## Label the correlation between Cog and SO
          Cognitive ~~ cor*Somatic
          ## Label the error variances of Per and SC
          Performance ~~ var_Per*Performance
          Self_confidence ~~ var_SC*Self_confidence"

## Display the model
plot(model, layout="spring")

RAM <- lavaan2RAM(model, obs.variables=c("Performance", "Cognitive",
                                         "Somatic", "Self_confidence"))
RAM

A1 <- RAM$A
S1 <- RAM$S
################################################################################

## Second stage analysis
fixed2 <- tssem2(fixed1, Amatrix=A1, Smatrix=S1, diag.constraints=TRUE,
                 intervals.type="LB", model.name="TSSEM2 Becker09",
                 mx.algebras=list( Cog=mxAlgebra(Cog2SC*SC2Per, name="Cog"),
                                   SO=mxAlgebra(SO2SC*SC2Per, name="SO"),
                                   Cog_SO=mxAlgebra(Cog2SC*SC2Per+SO2SC*SC2Per,
                                   name="Cog_SO")) )
summary(fixed2)

## Display the model with the parameter estimates
plot(fixed2, layout="spring")

#### Fixed-effects model: with type of sport as cluster
## First stage analysis
cluster1 <- tssem1(Becker09$data, Becker09$n, method="FEM",
                   cluster=Becker09$Type_of_sport)
summary(cluster1)

## Second stage analysis
cluster2 <- tssem2(cluster1, Amatrix=A1, Smatrix=S1, diag.constraints=TRUE,
                 intervals.type="LB", model.name="TSSEM2 Becker09",
                 mx.algebras=list( Cog=mxAlgebra(Cog2SC*SC2Per, name="Cog"),
                                   SO=mxAlgebra(SO2SC*SC2Per, name="SO"),
                                   Cog_SO=mxAlgebra(Cog2SC*SC2Per+SO2SC*SC2Per,
                                   name="Cog_SO")) )
summary(cluster2)

## Convert the model to semPlotModel object with 2 plots
## Use the short forms of the variable names
my.plots <- lapply(X=cluster2, FUN=meta2semPlot, manNames=c("Per","Cog","SO","SC") )

## Load the library
library("semPlot")

## Setup two plots
layout(t(1:2))
## The labels are overlapped. We may modify it by using layout="spring"
semPaths(my.plots[[1]], whatLabels="est", nCharNodes=10, color="orange",
         layout="spring", edge.label.cex=0.8)
title("Individual sport")
semPaths(my.plots[[2]], whatLabels="est", nCharNodes=10, color="skyblue",
         layout="spring", edge.label.cex=0.8)
title("Team sport")


#### Random-effects model
## First stage analysis
random1 <- tssem1(Becker09$data, Becker09$n, method="REM", RE.type="Diag")
summary(random1)

## Second stage analysis
random2 <- tssem2(random1, Amatrix=A1, Smatrix=S1, diag.constraints=TRUE,
                  intervals.type="LB", model.name="TSSEM2 Becker09",
                  mx.algebras=list( Cog=mxAlgebra(Cog2SC*SC2Per, name="Cog"),
                                    SO=mxAlgebra(SO2SC*SC2Per, name="SO"),
                                    Cog_SO=mxAlgebra(Cog2SC*SC2Per+SO2SC*SC2Per,
                                    name="Cog_SO")) )
summary(random2)

## Display the model
plot(random2, what="path", layout="spring")

## Display the model with the parameter estimates
plot(random2, layout="spring", color="yellow")

#### Univariate r approach
#### First stage of the analysis
uni1 <- uniR1(Becker09$data, Becker09$n)
uni1

#### Second stage of analysis using OpenMx
model2 <- "## Regression paths
           Performance ~ Cog2Per*Cognitive + SO2Per*Somatic + SC2Per*Self_confidence
           Self_confidence ~ Cog2SC*Cognitive + SO2SC*Somatic
           ## Provide starting values for Cog and SO
           Cognitive ~~ start(1)*Cognitive
           Somatic ~~ start(1)*Somatic
           ## Label the correlation between Cog and SO
           Cognitive ~~ cor*Somatic
           ## Label the error variances of Per and SC
           Performance ~~ var_Per*Performance
           Self_confidence ~~ var_SC*Self_confidence"

RAM2 <- lavaan2RAM(model2, obs.variables=c("Performance", "Cognitive",
                                           "Somatic", "Self_confidence"))
RAM2

uni2mx <- uniR2mx(uni1, RAM=RAM2)
summary(uni2mx)

#### Second stage of analysis Using lavaan
model3 <- "## Regression paths
           Performance ~ Cognitive + Somatic + Self_confidence
           Self_confidence ~ Cognitive + Somatic"

uni2lavaan <- uniR2lavaan(uni1, model3)
lavaan::summary(uni2lavaan)
}
}
\keyword{datasets}
