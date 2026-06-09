# =============================================================================
# Run the CAP Methods Explorer on your LOCAL machine.
#
# The Shiny app must run on the same machine as your browser. Every method (and
# its compiled C++ kernel) comes from the CAPmethods package (see R/pkg_methods.R),
# so this script just makes sure the UI packages and CAPmethods are installed,
# then launches. The kernels are compiled once, when CAPmethods is installed.
#
# Usage (from the cloned repository root):
#     Rscript app/run_local.R            # uses port 7700
#     Rscript app/run_local.R 8080       # custom port
# or, in an R session opened anywhere:
#     source(".../CAPmethods/app/run_local.R")
#
# Requirements: R; internet access the first time (to install packages); a C++
# toolchain only if CAPmethods must be built from source here (macOS: Xcode
# Command Line Tools — `xcode-select --install`; Windows: Rtools; Linux: g++).
# =============================================================================

## ---- locate paths ----------------------------------------------------------
.args_all <- commandArgs(trailingOnly = FALSE)
.f <- sub("^--file=", "", .args_all[grep("^--file=", .args_all)])
app_dir <- if (length(.f)) dirname(normalizePath(.f)) else {
  d <- file.path(getwd(), "app")
  if (file.exists(file.path(d, "app.R"))) d else getwd()
}
root <- normalizePath(file.path(app_dir, ".."))
.port <- { p <- commandArgs(trailingOnly = TRUE)
           if (length(p) && !is.na(as.integer(p[1]))) as.integer(p[1]) else 7700L }
message("App dir: ", app_dir, "\nProject root: ", root, "\nPort: ", .port)

## ---- 1. install missing UI / runtime packages ------------------------------
# CAPmethods (installed below) pulls in the method-side deps (MASS, glmnet, nlme,
# brglm2, foreach, doParallel, Rcpp). mvtnorm is used only by the built-in
# example generators, so it is listed explicitly here.
pkgs <- c("shiny", "bslib", "bsicons", "DT", "plotly", "shinycssloaders",
          "markdown", "mvtnorm")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  message("Installing: ", paste(miss, collapse = ", "))
  install.packages(miss, repos = "https://cloud.r-project.org")
}

## ---- 2. ensure the CAPmethods package is installed -------------------------
# The app gets every method from CAPmethods. If it is missing, install it from
# the repository source (the app folder's parent has the package DESCRIPTION);
# this builds the C++ kernels once. Falls back to GitHub if the source is absent.
if (!requireNamespace("CAPmethods", quietly = TRUE)) {
  if (file.exists(file.path(root, "DESCRIPTION"))) {
    message("Installing CAPmethods from the repository source (", root,
            ") — builds kernels once ...")
    install.packages(root, repos = NULL, type = "source")
  } else {
    message("Installing CAPmethods from GitHub ...")
    if (!requireNamespace("remotes", quietly = TRUE))
      install.packages("remotes", repos = "https://cloud.r-project.org")
    remotes::install_github("zhaoyi1026/CAPmethods")
  }
}
if (!requireNamespace("CAPmethods", quietly = TRUE))
  stop("CAPmethods could not be installed — the app needs it for every method.")
message("CAPmethods ", as.character(utils::packageVersion("CAPmethods")), " ready.")

## ---- 3. launch (opens your default browser) --------------------------------
message("Launching app — your browser should open at http://localhost:", .port)
shiny::runApp(app_dir, port = .port, launch.browser = TRUE, host = "127.0.0.1")
