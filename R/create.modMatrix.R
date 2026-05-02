#' Create a moderator matrix used in OSMASEM
#' 
#' It creates a moderator matrix used in OSMASEM.
#' 
#' 
#' @param RAM A RAM object including a list of matrices of the model returned
#' from \code{\link[metaSEM]{lavaan2RAM}}.
#' @param output Whether the output is an "A" or "S" matrix.
#' @param mod A string of moderator in the dataset.
#' @return A character matrix.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @keywords utilities
#' @examples
#' 
#' ## A multiple regression model
#' model <- "y ~ x1 + x2
#'           x1 ~~ 1*x1
#'           x2 ~~ 1*x2
#'           x1 ~~ x2"
#' 
#' ## RAM specification
#' RAM <- lavaan2RAM(model, obs.variables=c("y", "x1", "x2"))
#' 
#' ## Create a moderator matrix on A with "meanAge as the moderator.
#' A1 <-  create.modMatrix(RAM=RAM, output="A", mod="meanAge")
#' A1
#' 
#' ## Create a moderator matrix on S with "meanAge as the moderator.
#' S1 <-  create.modMatrix(RAM=RAM, output="S", mod="meanAge")
#' S1
#' 
create.modMatrix <- function(RAM, output=c("A", "S"), mod) {
    output <- match.arg(output)

    switch(output,
        A = { out <- RAM$A
        out[grep("\\*", out)] <- paste0("0*data.", mod)},
        S = { out <- RAM$S
        out[grep("\\*", out)] <- paste0("0*data.", mod)
        Diag(out) <- "0"})

    out
}

