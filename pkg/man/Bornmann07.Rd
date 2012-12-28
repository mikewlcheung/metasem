\name{Bornmann07}
\alias{Bornmann07}
\docType{data}
\title{A Dataset from Bornmann et al. (2007)}
\description{
  A dataset from Bornmann et al. (2007) for three-level meta-analysis.
}
\usage{data(Bornmann07)}
\format{
}
\details{
    The variables are:
  \describe{
  \item{ID}{ID of the study}
  \item{Study}{Study name}
  \item{Cluster}{Cluster for effect sizes}
  \item{logOR}{Effect size: log odds ratio}
  \item{v}{Sampling variance of logOR}
  \item{Year}{Year of publication}
  \item{Type}{Type of proposal: either \bold{Grants} or \bold{Fellowship}}
  \item{Discipline}{Discipline of proposal: either \bold{Physical sciences}, \bold{Life sciences/biology}, 
	\bold{Social sciences/humanities} or \bold{Multidisciplinary})}
  \item{Country}{Country of proposal: either \bold{United States},
	\bold{Canada}, \bold{Australia}, \bold{United Kingdom} or \bold{Europe}}
}
}
\source{
  Bornmann, L., Mutz, R., & Daniel, H.-D. (2007). Gender differences in grant peer review: A meta-analysis. \emph{Journal of Informetrics}, \bold{1(3)}, 226-238. doi:10.1016/j.joi.2007.03.001
}
\references{
  Marsh, H. W., Bornmann, L., Mutz, R., Daniel, H.-D., & O'Mara, A. (2009). Gender Effects in the Peer Reviews of Grant Proposals: A Comprehensive Meta-Analysis Comparing Traditional and Multilevel Approaches. \emph{Review of Educational Research}, \bold{79(3)}, 1290-1326. doi:10.3102/0034654309334143
}
\examples{
data(Bornmann07)

#### ML estimation method
## No predictor
summary( meta3(y=logOR, v=v, cluster=Cluster, data=Bornmann07) )

## Type as a predictor
## Grants: 0
## Fellowship: 1
summary( meta3(y=logOR, v=v, x=(as.numeric(Type)-1),
               cluster=Cluster, data=Bornmann07) ) 

## Centered Year as a predictor
summary( meta3(y=logOR, v=v, x=scale(Year, scale=FALSE),
               cluster=Cluster, data=Bornmann07) )

#### REML estimation method
## No predictor
summary( reml3(y=logOR, v=v, cluster=Cluster, data=Bornmann07) )

## Type as a predictor
## Grants: 0
## Fellowship: 1
summary( reml3(y=logOR, v=v, x=(as.numeric(Type)-1),
               cluster=Cluster, data=Bornmann07) ) 

## Centered Year as a predictor
summary( reml3(y=logOR, v=v, x=scale(Year, scale=FALSE),
               cluster=Cluster, data=Bornmann07) )

## Handling missing covariates with FIML
## MCAR
## Set seed for replication
set.seed(1000000)

## Copy Bornmann07 to my.df
my.df <- Bornmann07
## "Fellowship": 1; "Grant": 0
my.df$Type_MCAR <- ifelse(Bornmann07$Type=="Fellowship", yes=1, no=0)

## Create 17 out of 66 missingness with MCAR
my.df$Type_MCAR[sample(1:66, 17)] <- NA
summary(meta3X(y=logOR, v=v, cluster=Cluster, x2=Type_MCAR, data=my.df))

## MAR
Type_MAR <- ifelse(Bornmann07$Type=="Fellowship", yes=1, no=0)

## Create 27 out of 66 missingness with MAR for cases Year<1996
index_MAR <- ifelse(Bornmann07$Year<1996, yes=TRUE, no=FALSE)
Type_MAR[index_MAR] <- NA

## Include auxiliary variable
summary(meta3X(y=logOR, v=v, cluster=Cluster, x2=Type_MAR, av2=Year, data=my.df))
}
\keyword{datasets}
