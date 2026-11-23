# Package load logic.
#
# Each CAP-family method defines identically-named internal helpers
# (capReg, cov.ls, gamma.solve, obj.func, ...). To keep them from colliding in
# the single package namespace, every method's code is sourced into its OWN
# private environment at load time (mirroring the Shiny app's design). The
# private environments enclose the package namespace, so the method code resolves
# imported functions (ginv, lme, rnorm, ...). Every .Call() in the method code
# names the package DLL (PACKAGE = .<prefix>_vN_dll, set to "CAPmethods" when the
# code is assembled by tools/build_caprpkg.R), so the registered kernels are used
# even if a method's source folder kernel (same symbol names) is loaded in the same
# session; the kernels carry per-method symbol prefixes
# (hd_/lcap_/mlcap_/coc_/med_/hdcov_/cluster_) so they never clash in the DLL.

# Registry of per-method environments (populated in .onLoad).
.cap_envs <- new.env(parent = emptyenv())

# Method id -> source file under inst/method/.
.cap_method_files <- c(
  hdcap     = "hdcap.R",
  lcap      = "lcap.R",
  mcap      = "mcap.R",
  coc       = "coc.R",
  mediation = "mediation.R",
  hdcov     = "hdcov.R",
  cluster   = "cluster.R"
)

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)

  for (m in names(.cap_method_files))
    assign(m, new.env(parent = ns), envir = .cap_envs)

  # MCAP needs the bundled von Mises-Fisher sampler in its env before its code.
  sys.source(system.file("method", "rvmf_function.R", package = pkgname),
             envir = get("mcap", envir = .cap_envs), keep.source = FALSE)

  for (m in names(.cap_method_files)) {
    f <- system.file("method", .cap_method_files[[m]], package = pkgname)
    sys.source(f, envir = get(m, envir = .cap_envs), keep.source = FALSE)
  }
}
