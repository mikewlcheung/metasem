.onAttach <- function(lib, pkg){
  libMatrix <- installed.packages()
  packageStartupMessage("Loading metaSEM version ", libMatrix["metaSEM", "Version"], ".")
  packageStartupMessage("You may refer to the vignette for the examples.\n")
}
