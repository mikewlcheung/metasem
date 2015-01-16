\name{pattern.n}
\alias{pattern.n}
\title{Display the Accumulative Sample Sizes for the Covariance Matrix
}
\description{It displays the accumulative sample sizes for the
  covariance matrix.
}
\usage{
pattern.n(x, n)
}
\arguments{
  \item{x}{A list of square matrices}
  \item{n}{A vector of sample sizes.
}
}
\value{A square matrix of the accumulative sample sizes of the
  input matrices.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## Show the pattern of missing data
pattern.n(Hunter83$data, Hunter83$n)

#             Ability Knowledge Work sample Supervisor
# Ability        3815      3372        3281       3605
# Knowledge      3372      3532        2998       3322
# Work sample    3281      2998        3441       3231
# Supervisor     3605      3322        3231       3765
}
\keyword{utilities}

