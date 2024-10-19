#' @importFrom basilisk BasiliskEnvironment
APL_env <- BasiliskEnvironment(
    envname = "APL-env",
    pkgname = "APL",
    packages = c(
        "numpy==1.26.4",
        "pytorch==2.3.0",
        "libstdcxx-ng==14.2.0",
        "libprotobuf==4.25.3",
        "mkl==2023.2.0"
    )
)
