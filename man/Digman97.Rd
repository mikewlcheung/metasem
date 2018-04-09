\name{Digman97}
\alias{Digman97}
\docType{data}
\title{Factor Correlation Matrices of Big Five Model from Digman (1997)}
\description{
The data set includes fourteen studies of the factor correlation matrices of the Five-Factor Model of personality reported by Digman (1997).
}
\usage{data(Digman97)}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of 14 studies of correlation matrices. The
	variables are \emph{Agreeableness} (A), \emph{Conscientiousness}
	(C), \emph{Emotional Stability} (ES), \emph{Extraversion} (E) and \emph{Intellect} (I)}
  \item{n}{A vector of sample sizes}
  \item{cluster}{Types of participants of the studies}
}
}
\source{
Digman, J.M. (1997). Higher-order factors of the Big Five. \emph{Journal of Personality and Social Psychology}, \bold{73}, 1246-1256.
}
\references{
 Cheung, M. W.-L., & Chan, W. (2005). Classifying correlation matrices into relatively homogeneous subgroups: A cluster analytic approach. \emph{Educational and Psychological Measurement}, \bold{65}, 954-979.
}
\examples{
\dontrun{
Digman97

##### Fixed-effects TSSEM
fixed1 <- tssem1(Digman97$data, Digman97$n, method="FEM")
summary(fixed1)

## Factor covariance among latent factors
Phi <- matrix(c(1,"0.3*cor","0.3*cor",1), ncol=2, nrow=2)

## Error covariance matrix
Psi <- Diag(c("0.2*e1","0.2*e2","0.2*e3","0.2*e4","0.2*e5"))

## S matrix
S1 <- bdiagMat(list(Psi, Phi))

## This step is not necessary but it is useful for inspecting the model.
dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- c("A","C","ES","E","I","Alpha","Beta")

## Display S1
S1

## A matrix
Lambda <-
matrix(c(".3*Alpha_A",".3*Alpha_C",".3*Alpha_ES",rep(0,5),".3*Beta_E",".3*Beta_I"),
       ncol=2, nrow=5)
A1 <- rbind( cbind(matrix(0,ncol=5,nrow=5), Lambda),
             matrix(0, ncol=7, nrow=2) )

## This step is not necessary but it is useful for inspecting the model.
dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- c("A","C","ES","E","I","Alpha","Beta")

## Display A1
A1

## F matrix to select the observed variables
F1 <- create.Fmatrix(c(1,1,1,1,1,0,0), as.mxMatrix=FALSE)

## Display F1
F1

################################################################################
## Alternative model specification in lavaan model syntax
model <- "## Factor loadings
          Alpha=~A+C+ES
          Beta=~E+I
          ## Factor correlation
          Alpha~~Beta"

## Display the model
plot(model)

RAM <- lavaan2RAM(model, obs.variables=c("A","C","ES","E","I"),
                  A.notation="on", S.notation="with")
RAM

A1 <- RAM$A
S1 <- RAM$S
F1 <- RAM$F

################################################################################
fixed2 <- tssem2(fixed1, Amatrix=A1, Smatrix=S1, Fmatrix=F1,
                 model.name="TSSEM2 Digman97")
summary(fixed2)

## Display the model with the parameter estimates
plot(fixed2)

#### Fixed-effects TSSEM with several clusters
#### Create a variable for different samples
#### Younger participants: Children and Adolescents
#### Older participants: others
cluster <- ifelse(Digman97$cluster \%in\% c("Children","Adolescents"),
                  yes="Younger participants", no="Older participants")

#### Show the cluster
cluster

## Example of Fixed-effects TSSEM with several clusters
fixed1.cluster <- tssem1(Digman97$data, Digman97$n, method="FEM",
                         cluster=cluster)
summary(fixed1.cluster)

fixed2.cluster <- tssem2(fixed1.cluster, Amatrix=A1, Smatrix=S1, Fmatrix=F1)
#### Please note that the estimates for the younger participants are problematic.
summary(fixed2.cluster)

## Setup two plots
layout(t(1:2))

## Plot the first group
plot(fixed2.cluster[[1]])
title("Younger participants")

## Plot the second group
plot(fixed2.cluster[[2]])
title("Older participants")

#### Random-effects TSSEM with random effects on the diagonals
random1 <- tssem1(Digman97$data, Digman97$n, method="REM",
                  RE.type="Diag")
summary(random1)

random2 <- tssem2(random1, Amatrix=A1, Smatrix=S1, Fmatrix=F1)
summary(random2)

## Display the model with the parameter estimates
plot(random2, color="green")
}
}
\keyword{datasets}
