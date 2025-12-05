pkgname <- "jmleIRT"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "jmleIRT-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('jmleIRT')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("biasCorrection")
### * biasCorrection

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: biasCorrection
### Title: Analytical First-Order Bias Correction for Joint Maximum
###   Likelihood Estimators in the Rasch Model
### Aliases: biasCorrection

### ** Examples

## Not run: 
##D # Example data matrix with 2 persons and 4 items:
##D X <- matrix(c(1, 0, 1, NA, 1, 0, 1, 1), nrow = 2, byrow = TRUE)
##D theta <- c(0.5, -0.5)
##D beta <- c(-0.2, 0.1, 0.3, -0.1)
##D I <- ncol(X)
##D corrected <- biasCorrectionJMLE(theta, beta, X, I)
##D print(corrected$theta)
##D print(corrected$beta)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("biasCorrection", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("prox_algorithm")
### * prox_algorithm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: prox_algorithm
### Title: PROX Estimation for the Rasch Model
### Aliases: prox_algorithm

### ** Examples

## Not run: 
##D data <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10, ncol = 10)
##D result <- prox_algorithm(data)
##D print(result$b)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("prox_algorithm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
