# =============================================================================
# CAP Methods Explorer -- Shiny app entry point
# -----------------------------------------------------------------------------
# A website to demonstrate and run the Covariate Assisted Principal (CAP)
# family of statistical methods. Architecture:
#   R/           generic, method-agnostic machinery (registry + modules)
#   methods/<id> one self-contained plugin per method (registers itself)
#
# To add a method: create methods/<id>/<id>_method.R that calls
# register_method(...). It is auto-discovered and gets its own page.
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(bsicons)
  library(DT)
  library(plotly)
  library(shinycssloaders)
  library(markdown)
})

# Under shiny::runApp() the working directory is the app directory.
APP_DIR <- getwd()

# ---- load generic machinery (registry FIRST: defines %||%, register_method) -
# local = TRUE keeps all definitions in this environment, where APP_DIR
# lives (runApp() evaluates app.R in a sandbox, not the global environment).
source(file.path(APP_DIR, "R", "registry.R"), local = TRUE)
for (f in c("io_helpers.R", "ui_helpers.R", "mod_method.R", "pkg_methods.R")) {
  source(file.path(APP_DIR, "R", f), local = TRUE)
}

# ---- discover & load method plugins ----------------------------------------
method_files <- list.files(file.path(APP_DIR, "methods"),
                           pattern = "_method\\.R$",
                           recursive = TRUE, full.names = TRUE)
for (mf in method_files) {
  tryCatch(source(mf, local = TRUE),
           error = function(e) warning(sprintf("Failed to load %s: %s",
                                               mf, conditionMessage(e))))
}
METHODS <- list_methods()

# ---- roadmap: the full CAP family (implemented + planned) -------------------
CAP_FAMILY <- list(
  list(id = "cap",            name = "CAP",            desc = "Original covariate assisted principal regression."),
  list(id = "hdcap",          name = "HDCAP",          desc = "High-dimensional CAP with covariance shrinkage & sparsity."),
  list(id = "cap-hdcov",      name = "CAP-HDcov",      desc = "CAP with high-dimensional predictors."),
  list(id = "lcap-invar",     name = "LCAP",           desc = "Longitudinal CAP, time-invariant projection with random effects."),
  list(id = "lcap-var",       name = "MCAP",           desc = "Multilevel CAP, cluster-varying projection."),
  list(id = "cap-coc",        name = "CAP-CoC",        desc = "Covariance-on-covariance regression."),
  list(id = "cap-mediation",  name = "CAP-mediation",  desc = "Mediation analysis with a graph mediator."),
  list(id = "cap-clustering", name = "CAP-clustering", desc = "Covariance matrix clustering.")
)

# =============================================================================
# Home page
# =============================================================================
home_ui <- function() {
  cards <- lapply(CAP_FAMILY, function(m) {
    impl <- !is.null(get_method(m$id))
    badge <- if (impl) status_badge(get_method(m$id)$status %||% "ready")
             else status_badge("planned")
    bslib::card(
      class = if (impl) "h-100 border-primary" else "h-100",
      bslib::card_header(tags$b(m$name), " ", badge),
      bslib::card_body(
        tags$p(class = "small", m$desc),
        if (impl)
          tags$span(class = "small text-success",
                    bsicons::bs_icon("arrow-right-circle"),
                    " Open from the top navigation.")
        else
          tags$span(class = "small text-muted", "Coming soon.")
      )
    )
  })
  bslib::page_fluid(
    bslib::card(
      bslib::card_body(
        h2("CAP Methods Explorer"),
        tags$p(class = "lead",
               "Demonstrate and run the Covariate Assisted Principal (CAP) family of methods for covariance-matrix outcomes."),
        tags$p("Each method has its own page where you can read an overview, run it on a built-in simulated example, or upload your own data and visualize the results. Select a method from the navigation bar above to begin."),
        tags$hr(),
        tags$p(class = "text-muted small",
               sprintf("%d of %d methods currently available.",
                       length(METHODS), length(CAP_FAMILY)))
      )
    ),
    h4("Methods"),
    do.call(bslib::layout_column_wrap,
            c(list(width = 1/2), cards))
  )
}

# =============================================================================
# Assemble navbar: Home + one page per implemented method
# =============================================================================
# order the tabs by the canonical CAP_FAMILY order (CAP, HDCAP, CAP-HDcov, LCAP,
# MCAP, ...); any registered method not in CAP_FAMILY is appended after.
.fam_order <- vapply(CAP_FAMILY, function(m) m$id, character(1))
.method_ids <- names(METHODS)
.ordered_ids <- c(intersect(.fam_order, .method_ids),
                  setdiff(.method_ids, .fam_order))
method_navs <- unname(lapply(.ordered_ids, function(id) {
  spec <- METHODS[[id]]
  bslib::nav_panel(
    title = spec$name,
    methodUI(spec$id, spec)
  )
}))

ui <- do.call(
  bslib::page_navbar,
  c(
    list(
      title = tagList(bsicons::bs_icon("diagram-3"), " CAP Explorer"),
      id = "main_nav",
      theme = bslib::bs_theme(version = 5, bootswatch = "flatly",
                              primary = "#2c7fb8"),
      bslib::nav_panel(title = "Home", icon = bsicons::bs_icon("house"),
                       home_ui())
    ),
    method_navs,
    list(
      bslib::nav_spacer(),
      bslib::nav_item(tags$a(href = "https://github.com",
                             "About", class = "nav-link"))
    )
  )
)

# =============================================================================
# Server: wire up each method's module
# =============================================================================
server <- function(input, output, session) {
  for (spec in METHODS) {
    methodServer(spec$id, spec)
  }
}

shinyApp(ui, server)
