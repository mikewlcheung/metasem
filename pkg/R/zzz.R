## Removed to improve performance

## .onAttach <- function(lib, pkg){
##   libMatrix <- installed.packages()
##   packageStartupMessage("Loaded OpenMx version ", libMatrix["OpenMx", "Version"], ".")
##   packageStartupMessage("Loading metaSEM version ", libMatrix["metaSEM", "Version"], ".")
##   packageStartupMessage("You may refer to the vignette for the examples.\n")
## }

.onAttach <- function(lib, pkg){
        ## Use NPSOL as the default optimizer "CSOLNP"
        default <- "NPSOL"
        mxOption(NULL, "Default optimizer", default)

        if (default=="NPSOL") {
                packageStartupMessage("\"NPSOL\" is set as the default optimizer in OpenMx.")
                packageStartupMessage("You may change it to \"CSOLNP\" by calling:")
                packageStartupMessage("mxOption(NULL, \"Default optimizer\", \"CSOLNP\")\n")
        } else {
                packageStartupMessage("\"CSOLNP\" is set as the default optimizer in OpenMx.")
                packageStartupMessage("You may change it to \"NPSOL\" by calling:")
                packageStartupMessage("mxOption(NULL, \"Default optimizer\", \"NPSOL\")\n")
        }
}
