.onAttach <- function(lib, pkg){
  libMatrix <- installed.packages()
  packageStartupMessage("Loaded OpenMx version ", libMatrix["OpenMx", "Version"], ".")
  packageStartupMessage("Loading metaSEM version ", libMatrix["metaSEM", "Version"], ".")
  packageStartupMessage("You may refer to the vignette for the examples.\n")
}
