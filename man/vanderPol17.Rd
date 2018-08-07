\name{vanderPol17}
\alias{vanderPol17}
\docType{data}
\title{Dataset on the effectiveness of multidimensional family therapy in treating
  adolescents with multiple behavior problems
}
\description{
  This dataset includes 61 effect sizes from 19 manuscripts nested from
  8 studies reported by van der Pol et al. (2017). It studies the
  effectiveness of multidimensional family therapy in treating
  adolescents with multiple behavior problems.
}
\usage{data(vanderPol17)}

\details{
A list of data with the following structure:
\describe{
  \item{Number}{Number of the effect size.}
  \item{Study}{Authors of the studies.}
  \item{N}{Total sample size.}
  \item{N_target}{Sample size in the target group.}
  \item{N_control}{Sample size in the control group.}
  \item{Comparison_condition}{Either cognitive behavioral therapy
	(\code{CBT}), combined treatment (\code{CT}) or group therapy (\code{Group}).}
  \item{Study_ID}{Level-3 cluster.}
  \item{Age_mean}{Mean age of the participants.}
  \item{Fllow_up}{Follow-up duration (in months).}
  \item{Per_Males}{Percentage of males.}
  \item{Per_Minorities}{Percentage of minorities.}
  \item{Per_Conduct_disorder}{Percentage of participants with conduct disorder}
  \item{Per_Servere_cannabis_users}{Percentage of participants of
	servere cannabis use.}
  \item{Outcome_measure}{Either substance abuse, delinquency, externalizing and internalizing psychopathology, and
family functioning}
  \item{d}{Effect size in Cohen's d.}
  \item{v}{Sampling variance of d}
}
}

\source{
  van der Pol, T. M., Hoeve, M., Noom, M. J., Stams, G. J. J. M., Doreleijers, T. A. H., van Domburgh, L., & Vermeiren, R. R. J. M. (2017). Research Review: The effectiveness of multidimensional family therapy in treating adolescents with multiple behavior problems - a meta-analysis. \emph{Journal of Child Psychology and Psychiatry}, \bold{58}(5), 532-545. https://doi.org/10.1111/jcpp.12685
}

\examples{
data(vanderPol17)
}
\keyword{datasets}
