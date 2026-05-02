#' Check the correctness of the RAM formulation
#' 
#' It provides simple checks on the correctness of the RAM formulation.
#' 
#' 
#' @param Amatrix An asymmetric matrix in the RAM specification with
#' \code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be
#' converted into \code{\link[OpenMx]{MxMatrix-class}} by the
#' \code{as.mxMatrix} function.
#' @param Smatrix A symmetric matrix in the RAM specification with
#' \code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be
#' converted into \code{\link[OpenMx]{MxMatrix-class}} by the
#' \code{as.mxMatrix} function.
#' @param cor.analysis Logical. Analysis of correlation or covariance
#' structure. There are additional checks for cor.analysis=\code{TRUE}.
#' @return It returns silently if no error has been detected; otherwise, it
#' returns a warning message.
#' @author Mike W.-L. Cheung <mikewlcheung@@nus.edu.sg>
#' @seealso \code{\link[metaSEM]{as.mxMatrix}},
#' \code{\link[metaSEM]{lavaan2RAM}}
#' @keywords utilities
#' @examples
#' 
#' \donttest{ 
#' ## Digman97 example
#' model1 <- "## Factor loadings
#'            Alpha=~A+C+ES
#'            Beta=~E+I
#'            ## Factor correlation
#'            Alpha~~Beta"
#' 
#' RAM1 <- lavaan2RAM(model1, obs.variables=c("A","C","ES","E","I"),
#'                    A.notation="on", S.notation="with")
#' RAM1
#' 
#' ## The model is okay.    
#' checkRAM(Amatrix=RAM1$A, Smatrix=RAM1$S)
#' 
#' ## Hunter83 example    
#' model2 <- "## Regression paths
#'            Job_knowledge ~ A2J*Ability
#'            Work_sample ~ A2W*Ability + J2W*Job_knowledge
#'            Supervisor ~ J2S*Job_knowledge + W2S*Work_sample
#' 
#'            ## Fix the variance of Ability at 1
#'            Ability ~~ 1*Ability
#' 
#'            ## Label the error variances of the dependent variables
#'            Job_knowledge ~~ VarE_J*Job_knowledge
#'            Work_sample ~~ VarE_W*Work_sample
#'            Supervisor ~~ VarE_S*Supervisor"
#' 
#' RAM2 <- lavaan2RAM(model2, obs.variables=c("Ability","Job_knowledge",
#'                    "Work_sample","Supervisor"))
#' 
#' ## The model is okay.     
#' checkRAM(Amatrix=RAM2$A, Smatrix=RAM2$S)   
#' }
#' 
checkRAM <- function(Amatrix, Smatrix, cor.analysis=TRUE) {
    if (missing(Amatrix)&missing(Smatrix)) {
        warning("Either 'Amatrix' or 'Smatrix' must be present.")
    }    

    ## Check A
    if (!missing(Amatrix)) {
        ## Convert A into mxMatrix if they are not yet.
        if (is.matrix(Amatrix)) {
            Amatrix <- as.mxMatrix(Amatrix, name="A")
        } else {
            Amatrix@name <- "A"
        }

        ## A_name <- deparse(substitute(A))

        ## Check diagonals: either free or non-zero values
        if ( any(Diag(Amatrix$free)) | any(Diag(Amatrix$values)!=0) ) {
            warning("Diagonals of the 'Amatrix' must be zeros.\n")
        }

        ## Check both lower tri & upper tri are TRUE
        if ( any(Amatrix$free & t(Amatrix$free)) ) {
            warning("Non-recursive models are not allowed in the 'Amatrix'.\n")
        }
    }

    ## Check S
    if (!missing(Smatrix)) {
        if (is.matrix(Smatrix)) {
            Smatrix <- as.mxMatrix(Smatrix, name="S")
        } else {
            Smatrix@name <- "S"
        }

        ## S_name <- deparse(substitute(S))
        ## Cannot check 'free' for definition variables
        ## Only check labels!!!
        S_labels <- Smatrix$labels

        ## Check symmetric
        if (!isSymmetric(Smatrix$free)) {
            warning("The free parameters of the 'Smatrix' must be symmetric.\n")
        }

        if (!isSymmetric(Smatrix$labels)) {
            warning("The labels of 'Smatrix' must be symmetric.\n")
        }

        if (!isSymmetric(Smatrix$values)) {
            warning("The values of 'Smatrix' must be symmetric.\n")
        }
        
        ## ## Check digaonals
        ## if (SdiagZero==TRUE & any(!is.na(Diag(S_labels)))) {
        ##     warning("Diagonal of 'S' must be 0.\n")
        ## }

        ## Check both A and S: Variances of IVs must be fixed at 1 and DVs must be free
        ## Limitation: it may still give warnings when there are DVs with fixed parameters in A
        if (cor.analysis==TRUE & !missing(Amatrix)) {
            ## Check A if it is a DV
            dv <- apply(Amatrix$free, 1, any)

            if ( any(diag(Smatrix$free)[!dv]) | !all(diag(Smatrix$values)[!dv]==1) ) {
                warning("The variances of the independent variables in 'Smatrix' must be fixed at 1.")
            }            
            if (!all(diag(Smatrix$free)[dv])) {
                warning("The variances of the dependent variables in 'Smatrix' should be free.")
            }                        
        }        
    }      
}        
