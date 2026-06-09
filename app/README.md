# CAP Methods Explorer (Shiny app)

A Shiny web app to **demonstrate** and **run** the Covariate-Assisted Principal
(CAP) family of methods for covariance-matrix outcomes. It is a companion to the
[**CAPmethods**](../) R package — every method (and its compiled C++ kernel) comes
from that package, so the app itself is pure R/UI code.

For each method the site offers:

- **Overview** — a plain-language explanation and a link to the paper.
- **Run / Demo** — run on a built-in simulated example *or* upload your own data,
  tune parameters, and get results.
- **Results** — interactive tables (estimates + inference) and plots (loadings,
  subject scores, …), all downloadable.

| Method | Notes |
|--------|-------|
| CAP | Covariate-assisted principal regression |
| HDCAP | High-dimensional CAP (covariance shrinkage) |
| CAP-HDcov | CAP with high-dimensional covariates (sparse + selection inference) |
| LCAP | Longitudinal CAP, time-invariant projection + random effects |
| MCAP | Multilevel CAP, cluster-varying projection (von Mises–Fisher) |
| CAP-CoC | Covariance-on-covariance regression |
| CAP-mediation | Mediation with a covariance mediator |
| CAP-clustering | Clustering of covariance matrices |

## Install and run locally

The app must run on the **same machine as your browser**.

**1. Get the code.** Clone (or download) the repository:

```bash
git clone https://github.com/zhaoyi1026/CAPmethods.git
cd CAPmethods
```

**2. Install the `CAPmethods` package** (the method engine — this compiles the C++
kernels once, so a C++ toolchain is needed only here, not to run the app):

```r
# from the repository root, install the package source you just cloned:
install.packages(".", repos = NULL, type = "source")
# ... or install the published version directly:
# remotes::install_github("zhaoyi1026/CAPmethods")
```

A C++ toolchain is required to build `CAPmethods` from source — macOS:
`xcode-select --install`; Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/);
Linux: `g++`.

**3. Install the app's UI packages:**

```r
install.packages(c("shiny", "bslib", "bsicons", "DT", "plotly",
                   "shinycssloaders", "markdown", "mvtnorm"))
```

**4. Launch the app.** Either run the launcher (which checks/install steps 2–3 for
you and opens your browser):

```bash
Rscript app/run_local.R          # add a port number to override 7700
```

or, from an R session at the repository root:

```r
shiny::runApp("app", launch.browser = TRUE)
```

The app opens at `http://localhost:7700` (or the port you chose).

> The C++ kernels live in the installed `CAPmethods` package, so the app builds
> and loads nothing itself. The compiled binary is OS / CPU / R-version specific —
> don't copy an installed library between machines; reinstall instead.

## Architecture

The site is a **plugin host**: a generic Shiny module renders whatever each method
*declares*, so methods are independent.

```
app/
├── app.R                  # entry point: loads machinery, discovers plugins,
│                          #   builds the navbar (Home + one page per method)
├── R/                     # generic, method-agnostic machinery
│   ├── registry.R         #   register_method() / list_methods() / get_method()
│   ├── io_helpers.R       #   CSV parsing: long Y -> list of matrices, X -> matrix
│   ├── ui_helpers.R       #   param-spec -> Shiny inputs; table styling
│   ├── mod_method.R       #   the generic Shiny module that renders ANY method
│   └── pkg_methods.R      #   cap_pkg_env(): pull a method's fns from CAPmethods
└── methods/               # one self-contained plugin per method
    └── cap/
        ├── cap_method.R   #   cap_pkg_env("hdcap") + register_method(...)
        ├── explain.md     #   the Overview text
        └── sample_data/   #   downloadable example CSVs
```

### Adding a method

Create `methods/<id>/<id>_method.R` and call `register_method()` with a spec; it is
auto-discovered (any `*_method.R`) and gets its own page. Use
`methods/cap/cap_method.R` as the template. Key spec fields: `id`/`name`/
`full_name` (identity), `data_inputs` + `params` (the upload/parameter UI),
`example()` (built-in data, optionally with `truth`), `parse(files, opts)`
(uploaded data → internal dataset), `run(d, params)` (run the method),
`summarize(res)` (result tables) and `plots(res)` (result plots).

Each plugin pulls its method's functions from the installed `CAPmethods` package
via `cap_pkg_env()`:

```r
# CAPmethods method ids: hdcap, lcap, mcap, coc, mediation, hdcov, cluster
.hdcap_env <- cap_pkg_env("hdcap")          # capReg(), cap_beta_boot(), ...
fit <- .hdcap_env$capReg(Y, X, nD = 1, cov.shrinkage = FALSE)
```
