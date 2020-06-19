\name{create.modMatrix}
\alias{create.modMatrix}
\title{Create a moderator matrix used in OSMASEM
}
\description{It creates a moderator matrix used in OSMASEM.  
}
\usage{
create.modMatrix(RAM, output=c("A", "S"), mod)
}
\arguments{
\item{RAM}{A RAM object including a list of matrices of the model
  returned from \code{\link[metaSEM]{lavaan2RAM}}.}
\item{output}{Whether the output is an "A" or "S" matrix.}
\item{mod}{A string of moderator in the dataset.}
}

\value{A character matrix.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\examples{
## A multiple regression model
model <- "y ~ x1 + x2
          x1 ~~ 1*x1
          x2 ~~ 1*x2
          x1 ~~ x2"

## RAM specification
RAM <- lavaan2RAM(model, obs.variables=c("y", "x1", "x2"))

## Create a moderator matrix on A with "meanAge as the moderator.
A1 <-  create.modMatrix(RAM=RAM, output="A", mod="meanAge")
A1

## Create a moderator matrix on S with "meanAge as the moderator.
S1 <-  create.modMatrix(RAM=RAM, output="S", mod="meanAge")
S1
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
