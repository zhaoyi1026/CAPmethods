# Assemble the CAPmethods R package under CAP-Rpkg/260605/CAPmethods
# - copies the 7 C++ kernels into src/ (unique filenames)
# - strips each method's kernel-loader/dependency block into inst/method/
# - generates src/init.c (native-routine registration)
#
# 0.2.0 (2026-11-23): sources = HDCAP/V9, LCAP_gamma-invar/V7, LCAP_gamma-var/V5,
# CAP-CoC/V4, CAP-mediation/V5, CAP-HDcov/V3, CAP-clustering/V5.  These name their
# kernel in every .Call (PACKAGE = .<prefix>_vN_dll, set by a local({...}) block
# that builds/dyn.loads the folder's kernel); here that block becomes
# .<prefix>_vN_dll <- "CAPmethods", so every .Call resolves to the routines
# registered in the package DLL.

ROOT <- "/Users/yz125/Dropbox/MyFolder/Biostat-IU/Projects/cap/code"
PKG  <- file.path(ROOT, "CAP-Rpkg", "260605", "CAPmethods")

dir.create(file.path(PKG, "R"),           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PKG, "src"),         recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PKG, "inst", "method"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PKG, "man"),         recursive = TRUE, showWarnings = FALSE)

## ---- 1. kernels -> src/ (unique names) -------------------------------------
kernels <- c(
  "HDCAP/V9/cap_kernels.cpp"            = "src/hd_kernels.cpp",
  "LCAP_gamma-invar/V7/cap_kernels.cpp"= "src/lcap_kernels.cpp",
  "LCAP_gamma-var/V5/mlcap_kernels.cpp"= "src/mlcap_kernels.cpp",
  "CAP-CoC/V4/coc_kernels.cpp"         = "src/coc_kernels.cpp",
  "CAP-mediation/V5/med_kernels.cpp"   = "src/med_kernels.cpp",
  "CAP-HDcov/V3/hdcov_kernels.cpp"     = "src/hdcov_kernels.cpp",
  "CAP-clustering/V5/cluster_kernels.cpp" = "src/cluster_kernels.cpp"
)
for (s in names(kernels))
  file.copy(file.path(ROOT, s), file.path(PKG, kernels[[s]]), overwrite = TRUE)

## ---- 2. strip the loader block out of each method file ---------------------
# Removes: dependency attach lines (library/require/ps/install), the
# ".<prefix>_dir <- function(){...}" helper, the "if(!is.loaded(...)){...}"
# kernel-load block, plus any extra brace-block whose first line matches a regex
# (used to drop MCAP's Rfast/rvmf attach block). Brace-matched so nested {} are
# handled and the real ".helper" functions that follow are preserved.
strip_loader <- function(lines, extra_block_re = character()) {
  remove <- logical(length(lines))

  dep_re <- paste0("^\\s*(library\\(|require\\(|",
                   "suppressWarnings\\(suppressMessages\\(try\\(library|",
                   "ps\\s*<-\\s*c\\(|lapply\\(ps|install\\.packages\\()")
  remove[grepl(dep_re, lines)] <- TRUE

  block_re <- c("^\\.[A-Za-z]+_v[0-9]+_dir\\s*<-\\s*function",  # .<prefix>_dir
                "if\\s*\\(\\s*!is\\.loaded",                     # kernel load
                "^\\s*if\\s*\\(\\s*!require",                    # if(!require()){install}
                extra_block_re)
  starts <- sort(unique(unlist(lapply(block_re, function(re) grep(re, lines)))))

  count <- function(s, ch) sum(gregexpr(ch, s, fixed = TRUE)[[1]] > 0)
  for (s in starts) {
    if (remove[s] && s %in% grep(dep_re, lines)) next  # already a single dep line
    depth <- 0L; started <- FALSE; i <- s
    repeat {
      depth <- depth + count(lines[i], "{") - count(lines[i], "}")
      remove[i] <- TRUE
      if (count(lines[i], "{") > 0) started <- TRUE
      if (started && depth <= 0L) break
      i <- i + 1L
      if (i > length(lines)) break
    }
  }
  lines[!remove]
}

# Replace each ".<prefix>_vN_dll <- local({...})" block (the folder version's
# kernel build + dyn.load, whose value names the loaded DLL) by
# ".<prefix>_vN_dll <- \"CAPmethods\"": in the package the kernels are the
# routines registered in the package DLL.  Brace-matched, as strip_loader.
replace_dll_block <- function(lines) {
  re <- "^(\\.[A-Za-z]+_v[0-9]+_dll)\\s*<-\\s*local\\(\\{"
  count <- function(s, ch) sum(gregexpr(ch, s, fixed = TRUE)[[1]] > 0)
  repeat {
    s <- grep(re, lines)[1]
    if (is.na(s)) break
    var <- sub(paste0(re, ".*$"), "\\1", lines[s])
    depth <- 0L; i <- s
    repeat {
      depth <- depth + count(lines[i], "{") - count(lines[i], "}")
      if (depth <= 0L) break
      i <- i + 1L
      if (i > length(lines)) stop("unbalanced braces after ", var)
    }
    new <- c(sprintf("%s <- \"CAPmethods\"   # CAPmethods: the kernels are registered in the package DLL", var))
    lines <- c(lines[seq_len(s - 1L)], new, lines[-seq_len(i)])
  }
  lines
}

methods <- list(
  list(src = "HDCAP/V9/CAP_HD.R",                       dst = "hdcap.R"),
  list(src = "LCAP_gamma-invar/V7/Longitudinal_HDCAP.R",dst = "lcap.R"),
  # MCAP: drop the sampler choice (Rfast::rvmf or the bundled file); the package
  # always sources the bundled rvmf_function.R into the mcap env (zzz.R)
  list(src = "LCAP_gamma-var/V5/LCAP_varGamma.R",       dst = "mcap.R",
       extra = "^if\\s*\\(identical\\(getOption\\(\"mlcap\\.rvmf\"\\)"),
  list(src = "CAP-CoC/V4/COCReg.R",                     dst = "coc.R"),
  list(src = "CAP-mediation/V5/CAPMediation.R",         dst = "mediation.R"),
  list(src = "CAP-HDcov/V3/HCAP_Code.R",                dst = "hdcov.R"),
  list(src = "CAP-clustering/V5/CAP_Cluster.R",         dst = "cluster.R")
)
for (m in methods) {
  L <- readLines(file.path(ROOT, m$src), warn = FALSE)
  out <- replace_dll_block(L)
  out <- strip_loader(out, if (!is.null(m$extra)) m$extra else character())
  left <- grep("dyn\\.load|build_kernels|getLoadedDLLs|_dir\\(\\)|is\\.loaded|requireNamespace\\(\"Rfast\"", out, value = TRUE)
  left <- left[!grepl("^\\s*#", left)]
  if (length(left)) stop(m$dst, ": loader code left after stripping:\n", paste(left, collapse = "\n"))
  writeLines(out, file.path(PKG, "inst", "method", m$dst))
  cat(sprintf("  %-12s %4d -> %4d lines\n", m$dst, length(L), length(out)))
}
# MCAP's vMF sampler, sourced into the mcap env before mcap.R
file.copy(file.path(ROOT, "LCAP_gamma-var/V5/rvmf_function.R"),
          file.path(PKG, "inst", "method", "rvmf_function.R"), overwrite = TRUE)

## ---- 3. generate src/init.c (native routine registration) ------------------
routines <- c(
  hd_covls_cpp=1, hd_score_cpp=2, hd_capbeta_cpp=8, hd_capd1_cpp=9, hd_objfunc_cpp=6,
  lcap_covls_cpp=1, lcap_gammasolve_cpp=2, lcap_score_cpp=2, lcap_recar_cpp=16,
  lcap_capd1_cpp=15, lcap_objfunc_cpp=12,
  mlcap_compute_scores_cpp=3, mlcap_obj_func_cpp=15,
  coc_cov_ls_cpp=1, coc_cov_sk_x_cpp=1, coc_cov_sk_y_cpp=3, coc_score_cpp=2,
  coc_accum_cpp=2, coc_eigen_solve_cpp=2,
  med_score_cpp=2, med_accum_cpp=2, med_eigen_solve_cpp=2, med_cov_cpp=1,
  hdcov_precompute_cpp=1, hdcov_score_cpp=2, hdcov_accum_cpp=2,
  cluster_smat_cpp=1, cluster_score_cpp=2, cluster_accum_cpp=2, cluster_eigen_solve_cpp=2
)
decl <- vapply(names(routines), function(nm) {
  sexps <- paste(rep("SEXP", routines[[nm]]), collapse = ", ")
  sprintf("extern SEXP %s(%s);", nm, sexps)
}, character(1))
entries <- vapply(names(routines), function(nm)
  sprintf('    {"%s", (DL_FUNC) &%s, %d},', nm, nm, routines[[nm]]),
  character(1))

init_c <- c(
  "// Auto-generated native-routine registration for the CAPmethods package.",
  "#include <R.h>",
  "#include <Rinternals.h>",
  "#include <stdlib.h>",
  "#include <R_ext/Rdynload.h>",
  "",
  decl,
  "",
  "static const R_CallMethodDef CallEntries[] = {",
  entries,
  "    {NULL, NULL, 0}",
  "};",
  "",
  "void R_init_CAPmethods(DllInfo *dll) {",
  "    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);",
  "    R_useDynamicSymbols(dll, FALSE);",
  "    R_forceSymbols(dll, FALSE);",
  "}"
)
writeLines(init_c, file.path(PKG, "src", "init.c"))
cat("Wrote src/init.c with", length(routines), "routines\n")
cat("DONE assembling", PKG, "\n")
