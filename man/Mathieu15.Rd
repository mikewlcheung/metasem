\name{Mathieu15}
\alias{Mathieu15}
\docType{data}
\title{Correlation Matrices from Mathieu et al. (2015)}
\description{
The data set includes a list of correlation matrices of panel studies between
cohesion (C) and performance (P) in Mathieu et al. (2015,  Table 1). 
}
\usage{data(Mathieu15)}

\details{
A list of data with the following structure:
\describe{
  \item{data}{A list of studies of correlation matrices. The
	variables are \emph{C1}, \emph{P1}, \emph{C2}, and \emph{P2}.}
  \item{n}{A vector of sample sizes.}
  \item{Year}{Year of publication.}
  \item{Sample}{Sample characteristics.}
  \item{Student}{Whether the samples are student or non-student based on
	\code{Sample}.}
}
}
\source{
  Mathieu, J. E., Kukenberger, M. R., D'Innocenzo, L., & Reilly, G. (2015). Modeling reciprocal team cohesion-performance relationships, as impacted by shared leadership and members' competence. \emph{Journal of Applied Psychology}, \bold{100}(3), 713-734. https://doi.org/10.1037/a0038898
}
\examples{
\dontrun{
# TSSEM
## Model 1: no constraint
## Stage 1 analysis
tssem1.fit <- tssem1(Mathieu15$data, Mathieu15$n)
summary(tssem1.fit)

## Proposed model in lavaan syntax
model1 <- 'C2 ~ c2c*C1 + p2c*P1
           P2 ~ c2p*C1 + p2p*P1
           C1 ~~ c1withp1*P1
           C1 ~~ 1*C1
           P1 ~~ 1*P1
           C2 ~~ c2withp2*P2'

## Convert the lavaan model to RAM specification
RAM1 <- lavaan2RAM(model1, obs.variables=c("C1", "P1", "C2", "P2"))
RAM1

## Stage 2 analysis
tssem1b.fit <- tssem2(tssem1.fit, RAM=RAM1)
summary(tssem1b.fit)

plot(tssem1b.fit, col="yellow", edge.label.position=0.58)

## Model 2: Equality constraints on the path coefficient
## Proposed model with equal effects time 1 to time 2
model2 <- 'C2 ~ same*C1 + diff*P1
           P2 ~ diff*C1 + same*P1
           C1 ~~ c1withp1*P1
           C1 ~~ 1*C1
           P1 ~~ 1*P1
           C2 ~~ c2withp2*P2'

## Convert the lavaan model to RAM specification
RAM2 <- lavaan2RAM(model2, obs.variables=c("C1", "P1", "C2", "P2"))
RAM2
    
## Stage 2 analysis
tssem2b.fit <- tssem2(tssem1.fit, RAM=RAM2)
summary(tssem2b.fit)

## Compare the models with and without the constraints. 
anova(tssem1b.fit, tssem2b.fit)

## Plot the model
plot(tssem2b.fit, col="yellow", edge.label.position=0.60)


## OSMASEM
my.df <- Cor2DataFrame(Mathieu15)
    
head(my.df$data)

## Model without any moderator
osmasem.fit1 <- osmasem(model.name="No moderator", RAM=RAM1, data=my.df)
summary(osmasem.fit1)

## Extract the heterogeneity variance-covariance matrix
diag(VarCorr(osmasem.fit1))

plot(osmasem.fit1, col="yellow", edge.label.position=0.6)

## Model with student sample as a moderator on the regression coefficients
A1 <- create.modMatrix(RAM1, output="A", "Student")
A1

## Model with a moderator    
osmasem.fit2 <- osmasem(model.name="Student sample as a moderator", RAM=RAM1, 
                        Ax=A1, data=my.df)
summary(osmasem.fit2)

## Compare the models with and without the moderator
anova(osmasem.fit2, osmasem.fit1)

## Get the R2 of the moderator
osmasemR2(osmasem.fit2, osmasem.fit1)
}
}
\keyword{datasets}
