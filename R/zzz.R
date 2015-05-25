## Removed to improve performance

## .onAttach <- function(lib, pkg){
##   libMatrix <- installed.packages()
##   packageStartupMessage("Loaded OpenMx version ", libMatrix["OpenMx", "Version"], ".")
##   packageStartupMessage("Loading metaSEM version ", libMatrix["metaSEM", "Version"], ".")
##   packageStartupMessage("You may refer to the vignette for the examples.\n")
## }

.onAttach <- function(lib, pkg){
        ## SLSQP (not NPSOL) is used as the default optimizer.
        ## "central" (not "forward") is used as the default "Gradient algorithm".
        mxOption(NULL, "Gradient algorithm", "central")

        packageStartupMessage("\"SLSQP\" is set as the default optimizer in OpenMx.")
        packageStartupMessage("If \"SLSQP\" does not work well for you, e.g., there are many error codes,")
        packageStartupMessage("you may install the \"NPSOL\" optimizer from the OpenMx website and use it by calling:")
        packageStartupMessage("mxOption(NULL, \"Default optimizer\", \"NPSOL\")\n")
}
