\name{pattern.na}
\alias{pattern.na}
\title{Display the Pattern of Missing Data of a List of Square Matrices
}
\description{It displays the pattern of missing data (or pattern of data
  that are present) of a list of square matrices with the same dimensions.
}
\usage{
pattern.na(x, show.na = TRUE)
}
\arguments{
  \item{x}{A list of square matrices}
  \item{show.na}{If it is \code{TRUE}, it shows the pattern of missing
	data. If it is \code{FALSE}, it shows the pattern of data that are present.
}
}
\value{A square matrix of numerical values with the same dimensions of the
  input matrices.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## Show the pattern of missing data
pattern.na(Hunter83$data, show.na=TRUE)

#             Ability Knowledge Work sample Supervisor
# Ability           1         3           3          2
# Knowledge         3         2           4          3
# Work sample       3         4           2          3
# Supervisor        2         3           3          1

## Show the pattern of data that are present
pattern.na(Hunter83$data, show.na=FALSE)

#             Ability Knowledge Work sample Supervisor
# Ability          13        11          11         12
# Knowledge        11        12          10         11
# Work sample      11        10          12         11
# Supervisor       12        11          11         13
}
\keyword{utilities}

