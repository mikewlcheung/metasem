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

    mxOption(NULL, "Default optimizer", "SLSQP")
    mxOption(NULL, "Gradient algorithm", "central")
    mxOption(NULL, "Optimality tolerance", "6.3e-14")
    mxOption(NULL, "Gradient iterations", 2)

    packageStartupMessage('"SLSQP" is set as the default optimizer in OpenMx.')
    packageStartupMessage('mxOption(NULL, "Gradient algorithm") is set at "', mxOption(NULL, "Gradient algorithm"), '".')
    packageStartupMessage('mxOption(NULL, "Optimality tolerance") is set at "', mxOption(NULL, "Optimality tolerance"), '".')
    packageStartupMessage('mxOption(NULL, "Gradient iterations") is set at "', mxOption(NULL, "Gradient iterations"), '".')

    #packageStartupMessage('If "SLSQP" does not work well for you, e.g., there are many error codes,')
    #packageStartupMessage('you may install the "NPSOL" optimizer from the OpenMx website and use it by calling:')
    #packageStartupMessage('mxOption(NULL, "Default optimizer", "NPSOL")')
}

## Added .onLoad See https://github.com/mikewlcheung/metasem/issues/3 and https://github.com/OpenMx/OpenMx/issues/98
.onLoad <- function(lib, pkg){
    mxOption(NULL, "Default optimizer", "SLSQP")
    mxOption(NULL, "Gradient algorithm", "central")
    mxOption(NULL, "Optimality tolerance", "6.3e-14")
    mxOption(NULL, "Gradient iterations", 2)

##  packageStartupMessage('"SLSQP" is set as the default optimizer in OpenMx.')
##  packageStartupMessage('mxOption(NULL, "Gradient algorithm") is set at "', mxOption(NULL, "Gradient algorithm"), '".')
##  packageStartupMessage('mxOption(NULL, "Optimality tolerance") is set at "', mxOption(NULL, "Optimality tolerance"), '".')
##  packageStartupMessage('mxOption(NULL, "Gradient iterations") is set at "', mxOption(NULL, "Gradient iterations"), '".')
}
