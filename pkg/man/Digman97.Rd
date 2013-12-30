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
Digman97

##### Fixed-effects TSSEM
fixed1 <- tssem1(Digman97$data, Digman97$n, method="FEM")
summary(fixed1)

## S matrix
Phi <- matrix(c(1,"0.3*cor","0.3*cor",1), ncol=2, nrow=2)
S1 <- bdiagMat(list(Diag(c("0.2*e1","0.2*e2","0.2*e3","0.2*e4","0.2*e5")), Phi))

## This step is not necessary but it is useful for inspecting the model.
## dimnames(S1)[[1]] <- dimnames(S1)[[2]] <- c("A","C","ES","E","I","Alpha","Beta") 
## S1

S1 <- as.mxMatrix(S1)

## A matrix
Lambda <-
matrix(c(".3*Alpha_A",".3*Alpha_C",".3*Alpha_ES",rep(0,5),".3*Beta_E",".3*Beta_I"),
       ncol=2, nrow=5)
A1 <- rbind( cbind(matrix(0,ncol=5,nrow=5), Lambda),
             matrix(0, ncol=7, nrow=2) )

## This step is not necessary but it is useful for inspecting the model.
## dimnames(A1)[[1]] <- dimnames(A1)[[2]] <- c("A","C","ES","E","I","Alpha","Beta") 
## A1

A1 <- as.mxMatrix(A1)

## F matrix to select the observed variables
F1 <- create.Fmatrix(c(1,1,1,1,1,0,0), name="F1")

fixed2 <- tssem2(fixed1, Amatrix=A1, Smatrix=S1, Fmatrix=F1,
                 model.name="TSSEM2 Digman97")
summary(fixed2)

#### Fixed-effects TSSEM with several clusters
#### Create a variable for different samples
#### Younger participants: Children and Adolescents
#### Older participants: others
sample <- ifelse(Digman97$cluster \%in\% c("Children","Adolescents"),
                 yes="Younger participants", no="Older participants")

#### Show the sample
sample

## Example of Fixed-effects TSSEM with several clusters
fixed1.cluster <- tssem1(Digman97$data, Digman97$n, method="FEM", 
                         cluster=sample)
summary(fixed1.cluster)

fixed2.cluster <- tssem2(fixed1.cluster, Amatrix=A1, Smatrix=S1, Fmatrix=F1)
#### Please note that the estimates for the younger participants are problematic.
summary(fixed2.cluster)

#### Random-effects TSSEM with random effects on the diagonals
random1 <- tssem1(Digman97$data, Digman97$n, method="REM",
                  RE.type="Diag", I2="I2q")
summary(random1)

random2 <- tssem2(random1, Amatrix=A1, Smatrix=S1, Fmatrix=F1)
summary(random2)
}
\keyword{datasets}
