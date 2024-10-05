
#' @importFrom basilisk BasiliskEnvironment
APL_env <- BasiliskEnvironment(envname="APL-env",
    pkgname="APL",
    packages=c("numpy==1.26.4", "torch==2.2.0+cu118")
)
