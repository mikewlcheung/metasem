\name{osmasem}
\alias{osmasem}
\title{One-stage meta-analytic structural equation modeling
}
\description{It fits MASEM with the one-stage SEM (OSMASEM) approach.  
}
\usage{
osmasem(model.name="osmasem", Mmatrix, Tmatrix, data,
        intervals.type=c("z", "LB"), mxModel.Args=NULL,
        mxRun.Args=NULL, suppressWarnings=TRUE,
        silent=FALSE, run=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model.name}{A string for the model name in \code{\link[OpenMx]{mxModel}}. 
}
\item{Mmatrix}{A list of matrices of the model implied correlation matrix created by \code{create.vechsCor}.}
\item{Tmatrix}{A list of matrices of the heterogeneity variance-covariance matrix created by \code{create.Tau2}.}
\item{data}{A list of data created by \code{Cor2DataFrame}.}
\item{intervals.type}{Either \code{z} (default if missing) or
  \code{LB}. If it is \code{z}, it calculates the standard error. If it
  is \code{LB}, it calculates the 95\% likelihood-based CIs on the parameter estimates.}
\item{mxModel.Args}{A list of arguments which are passed to \code{\link[OpenMx]{mxModel}}.}
\item{mxRun.Args}{A list of arguments which are passed to \code{\link[OpenMx]{mxRun}}.}
  \item{suppressWarnings}{Logical. If \code{TRUE}, warnings are suppressed. The argument to be passed to \code{\link[OpenMx]{mxRun}}.}
  \item{silent}{Logical. The argument to be passed to \code{\link[OpenMx]{mxRun}}}
  \item{run}{Logical. If \code{FALSE}, only return the mx model without running the analysis.}   
  \item{\dots}{Not used yet.}
  
}

\value{An object of class \code{osmasem}}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{Cor2DataFrame}}, \code{\link[metaSEM]{create.vechsR}},
  \code{\link[metaSEM]{create.Tau2}}, \code{\link[metaSEM]{create.V}}
}

\examples{
\dontrun{
## Calculate the sampling covariance matrix of the correlation matrices.    
my.df <- Cor2DataFrame(Nohe15A1$data, Nohe15A1$n)

## Add the centered Lag as a moderator
my.df$data <- data.frame(my.df$data, Lag=scale(Nohe15A1$Lag, scale=FALSE),
                         check.names=FALSE)
head(my.df$data)

## Proposed model
model1 <- 'W2 ~ w2w*W1 + s2w*S1
           S2 ~ w2s*W1 + s2s*S1
           W1 ~~ w1WITHs1*S1
           W2 ~~ w2WITHs2*S2
           W1 ~~ 1*W1
           S1 ~~ 1*S1
           W2 ~~ Errw2*W2
           S2 ~~ Errs2*S2'
plot(model1)     

## Convert it into RAM specification    
RAM1 <- lavaan2RAM(model1, obs.variables=c("W1", "S1", "W2", "S2"))
RAM1

## Create vechs of the model implied correlation matrix
## with implicit diagonal constraints
M0 <- create.vechsCor(A0=RAM1$A, S0=RAM1$S, F0=NULL)

## Create heterogeneity variances
## RE.type= either "Diag" or "Symm"
## Transform= either "expLog" or "sqSD" for better estimation on variances
T0 <- create.Tau2(RAM=RAM1, RE.type="Diag")

## Fit the model
fit0 <- osmasem(model.name="No moderator", Mmatrix=M0, Tmatrix=T0, data=my.df)
summary(fit0)
VarCorr(fit0)
    
## Moderators with lag as a moderator on A matrix
A1 <- matrix(c(0,0,0,0,
               0,0,0,0,
               "0*data.Lag","0*data.Lag",0,0,
               "0*data.Lag","0*data.Lag",0,0),
             nrow=4, ncol=4, byrow=TRUE)
              
M1 <- create.vechsCor(A0=RAM1$A, S0=RAM1$S, F0=NULL,
                      Ax=A1, Sx=NULL)

fit1 <- osmasem(model.name="Lag as a moderator for Amatrix", Mmatrix=M1,
                Tmatrix=T0, data= my.df)
summary(fit1)
VarCorr(fit1)
    
anova(fit1, fit0)
osmasemR2(fit1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ osmasem }
