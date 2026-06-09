# =============================================================================
# Generic method module
# -----------------------------------------------------------------------------
# A single Shiny module that renders the full page for ANY registered method:
#   - Overview tab        : markdown explanation + manuscript link
#   - Run / Demo tab      : data source (example or upload) + parameters + run
#   - Results             : standardized tables + plots from the method spec
#
# Because everything is driven by the method `spec`, the same module powers all
# eight CAP-family methods.
# =============================================================================

methodUI <- function(id, spec) {
  ns <- NS(id)

  # ---- data source controls (example vs upload) ----------------------------
  upload_controls <- lapply(spec$data_inputs, function(d) {
    fileInput(ns(paste0("file_", d$id)), d$label,
              accept = d$accept %||%
                c(".rds", ".RData", ".rda", ".csv", ".tsv", ".txt"))
  })
  upload_help <- lapply(spec$data_inputs, function(d) {
    if (!is.null(d$help)) tags$li(tags$b(d$label), ": ", d$help)
  })

  sidebar <- bslib::sidebar(
    width = 340,
    radioButtons(ns("data_source"), "Data source",
                 choices = c("Built-in example" = "example",
                             "Upload my data"   = "upload"),
                 selected = "example"),
    conditionalPanel(
      sprintf("input['%s'] == 'upload'", ns("data_source")),
      tagList(
        upload_controls,
        if (isTRUE(spec$x_intercept_option))
          checkboxInput(ns("add_intercept"),
                        "Add intercept column to covariates", value = TRUE),
        tags$details(
          tags$summary("Expected file format"),
          tags$ul(class = "small", upload_help)
        )
      )
    ),
    conditionalPanel(
      sprintf("input['%s'] == 'example'", ns("data_source")),
      div(class = "alert alert-info p-2 small mb-2",
          spec$example_note %||%
            "A simulated dataset matching this method's model. Run to recover the truth."),
      if (is.function(spec$export_example) && length(spec$data_inputs) > 0)
        tagList(
          tags$small(class = "text-muted",
                     "Download this example to see the expected upload format:"),
          tags$div(lapply(spec$data_inputs, function(d)
            downloadButton(ns(paste0("dl_", d$id)),
                           paste("Download", d$label),
                           class = "btn-sm btn-outline-secondary mb-1 d-block")))
        )
    ),
    tags$hr(),
    tags$h6("Parameters"),
    make_param_panel(ns, spec$params),
    tags$hr(),
    actionButton(ns("run"), "Run analysis", class = "btn-primary w-100",
                 icon = icon("play"))
  )

  bslib::page_fluid(
    bslib::navset_card_tab(
      id = ns("tabs"),
      # ---------------------------------------------------------------- Overview
      bslib::nav_panel(
        title = "Overview", icon = bsicons::bs_icon("info-circle"),
        bslib::card(
          bslib::card_body(
            uiOutput(ns("overview"))
          )
        )
      ),
      # --------------------------------------------------------------- Run / Demo
      bslib::nav_panel(
        title = "Run / Demo", icon = bsicons::bs_icon("play-circle"),
        bslib::layout_sidebar(
          sidebar = sidebar,
          uiOutput(ns("data_preview")),
          shinycssloaders::withSpinner(
            uiOutput(ns("results_ui")), type = 6, color = "#2c7fb8"
          )
        )
      )
    )
  )
}

methodServer <- function(id, spec) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- Overview tab --------------------------------------------------------
    output$overview <- renderUI({
      md <- if (!is.null(spec$explain) && file.exists(spec$explain)) {
        HTML(markdown::markdownToHTML(spec$explain, fragment.only = TRUE))
      } else {
        tags$p(spec$short %||% "")
      }
      paper_link <- NULL
      if (!is.null(spec$paper) && !is.null(spec$paper$url)) {
        paper_link <- div(
          class = "alert alert-light border d-flex align-items-start gap-2",
          bsicons::bs_icon("journal-text", size = "1.4em"),
          div(
            tags$b("Reference: "),
            spec$paper$citation %||% spec$paper$title,
            tags$br(),
            tags$a(href = spec$paper$url, target = "_blank",
                   rel = "noopener noreferrer",
                   bsicons::bs_icon("box-arrow-up-right"), " ", spec$paper$url)
          )
        )
      }
      # withMathJax() loads MathJax and typesets the \( \) / $$ $$ delimiters
      # that markdownToHTML() emits.
      withMathJax(tagList(
        h3(spec$full_name), status_badge(spec$status),
        tags$hr(), md, paper_link
      ))
    })

    # ---- Downloadable example templates --------------------------------------
    if (is.function(spec$export_example)) {
      for (d in spec$data_inputs) {
        local({
          did <- d$id
          # data.frames download as CSV; native R objects (e.g. the response
          # list of matrices) download as .rds.
          output[[paste0("dl_", did)]] <- downloadHandler(
            filename = function() {
              obj <- spec$export_example(spec$example())[[did]]
              ext <- if (is.data.frame(obj)) "csv" else "rds"
              paste0(spec$id, "_", did, "_example.", ext)
            },
            content = function(file) {
              obj <- spec$export_example(spec$example())[[did]]
              if (is.data.frame(obj)) utils::write.csv(obj, file, row.names = FALSE)
              else saveRDS(obj, file)
            }
          )
        })
      }
    }

    # ---- Resolve the active dataset (example or uploaded) --------------------
    dataset <- reactive({
      if (input$data_source == "example") {
        validate(need(is.function(spec$example),
                      "No built-in example is available for this method."))
        spec$example()
      } else {
        # gather uploaded files and hand them to the method's parser
        files <- list()
        for (d in spec$data_inputs) {
          fi <- input[[paste0("file_", d$id)]]
          validate(need(!is.null(fi), paste0("Please upload: ", d$label)))
          files[[d$id]] <- read_upload(fi$datapath, fi$name)
        }
        opts <- list(add_intercept = if (isTRUE(spec$x_intercept_option))
          isTRUE(input$add_intercept) else TRUE)
        validate(need(is.function(spec$parse),
                      "This method cannot yet parse uploaded data."))
        spec$parse(files, opts)
      }
    })

    # ---- Data preview --------------------------------------------------------
    output$data_preview <- renderUI({
      d <- tryCatch(dataset(), error = function(e) {
        return(div(class = "alert alert-danger", conditionMessage(e)))
      })
      if (inherits(d, "shiny.tag")) return(d)
      bslib::card(
        class = "mb-3",
        bslib::card_header(bsicons::bs_icon("table"), " Data summary"),
        bslib::card_body(
          if (!is.null(d$preview_ui)) d$preview_ui else
            tags$p(spec$describe_data(d))
        )
      )
    })

    # ---- Run the analysis ----------------------------------------------------
    result <- eventReactive(input$run, {
      d <- dataset()
      params <- collect_params(input, spec$params)
      withProgress(message = "Running analysis...", value = 0.4, {
        spec$run(d, params)
      })
    })

    # ---- Results UI (tables + plots) -----------------------------------------
    output$results_ui <- renderUI({
      req(input$run)
      res <- tryCatch(result(), error = function(e) e)
      if (inherits(res, "error")) {
        return(div(class = "alert alert-danger",
                   tags$b("The analysis failed: "), conditionMessage(res)))
      }

      tables <- if (is.function(spec$summarize)) spec$summarize(res) else list()
      plots  <- if (is.function(spec$plots)) spec$plots(res) else list()

      # a right-aligned "Download CSV" button bound to a download output id
      dl_bar <- function(dlid) div(
        class = "d-flex justify-content-end mb-2",
        downloadButton(ns(dlid), "Download CSV",
                       class = "btn-sm btn-outline-secondary")
      )
      # register a CSV download handler for a captured data.frame
      register_csv <- function(dlid, df, label) local({
        d_ <- df; lab_ <- label
        output[[dlid]] <- downloadHandler(
          filename = function() paste0(spec$id, "_", safe_filename(lab_), ".csv"),
          content = function(file) utils::write.csv(d_, file, row.names = FALSE)
        )
      })

      # register dynamic table/plot outputs (each table downloadable as CSV)
      table_panels <- lapply(seq_along(tables), function(i) {
        nm <- names(tables)[i]
        outid <- paste0("tbl_", i); dlid <- paste0("dltbl_", i)
        output[[outid]] <- DT::renderDataTable(nice_table(tables[[i]]))
        register_csv(dlid, tables[[i]], nm)
        bslib::nav_panel(nm, dl_bar(dlid), DT::dataTableOutput(ns(outid)))
      })
      plot_panels <- lapply(seq_along(plots), function(i) {
        nm <- names(plots)[i]
        pl <- plots[[i]]
        outid <- paste0("plt_", i)
        # plots may carry a data.frame in pl$data -> offer it as CSV
        dl_ui <- NULL
        if (is.data.frame(pl$data)) {
          dlid <- paste0("dlplt_", i)
          register_csv(dlid, pl$data, nm)
          dl_ui <- dl_bar(dlid)
        }
        if (inherits(pl$plot, "plotly")) {
          output[[outid]] <- plotly::renderPlotly(pl$plot)
          bslib::nav_panel(nm, dl_ui, plotly::plotlyOutput(ns(outid), height = "460px"))
        } else {
          output[[outid]] <- renderPlot(pl$plot)
          bslib::nav_panel(nm, dl_ui, plotOutput(ns(outid), height = "460px"))
        }
      })

      do.call(bslib::navset_card_tab,
              c(list(title = tagList(bsicons::bs_icon("bar-chart-line"), " Results")),
                plot_panels, table_panels))
    })
  })
}
