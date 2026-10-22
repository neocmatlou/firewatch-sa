# =============================================================================
# FIREWATCH SA
# An Interactive Spatial Statistical Analysis of Satellite-Detected Fire
# Patterns in South Africa (2024-2025)
#
# Files in R/ (constants.R, analysis.R) are loaded automatically by Shiny.
# =============================================================================

library(shiny)
library(bslib)
library(bsicons)
library(leaflet)
library(sf)
library(terra)
library(spatstat.geom)
library(spatstat.explore)
library(ggplot2)

# ---- data -------------------------------------------------------------------
if (!file.exists("data/fires.rds"))
  stop("data/fires.rds not found - run scripts/01_prepare_data.R (with TEST_MODE <- FALSE) first.")

fires     <- readRDS("data/fires.rds")
provinces <- readRDS("data/provinces.rds")
windows   <- readRDS("data/windows.rds")
meta      <- readRDS("data/meta.rds")

year_choices   <- c(setNames(as.character(meta$years), meta$years), "Both years" = "all")
month_choices  <- c("All months" = "all", setNames(as.character(1:12), month.name))
method_choices <- c("Kernel density (KDE)" = "kde",
                    "Nearest neighbour (Clark-Evans)" = "nn",
                    "Ripley's K-function" = "ripley")

filter_fires <- function(region, year, month) {
  d <- fires
  if (region != "South Africa") d <- d[d$province == region, ]
  if (year  != "all") d <- d[d$year  == as.integer(year), ]
  if (month != "all") d <- d[d$month == as.integer(month), ]
  d
}

# ---- theme ------------------------------------------------------------------
fire_theme <- bs_theme(
  version = 5,
  bg = FIRE$ash, fg = "#231F1C",
  primary = FIRE$ember, secondary = FIRE$smoke,
  warning = FIRE$flame, danger = FIRE$deep,
  base_font    = font_google("Inter", local = FALSE),
  heading_font = font_google("Barlow Semi Condensed", wght = c(500, 600), local = FALSE),
  "border-radius" = "0.6rem"
)

# ---- UI ---------------------------------------------------------------------
ui <- page_navbar(
  title = tags$span(class = "brand", flame_svg(), "FIREWATCH", tags$span(class = "brand-sa", "SA")),
  window_title = "FIREWATCH SA",
  theme = fire_theme,
  fillable = "Explore",
  header = tagList(
    tags$link(rel = "stylesheet", href = "styles.css"),
    tags$div(class = "ember-line")
  ),

  # ---------------- Explore ----------------
  nav_panel(
    "Explore", value = "Explore", icon = bs_icon("map"),
    layout_sidebar(
      sidebar = sidebar(
        width = 310, title = "Analysis settings",
        selectInput("region", "Region", REGIONS),
        selectInput("year", "Year", year_choices, selected = "all"),
        selectInput("month", "Month", month_choices),
        radioButtons("method", "Spatial method", method_choices),

        conditionalPanel(
          "input.method == 'kde'",
          radioButtons("bw", "Bandwidth",
                       c("Scott's rule (data-driven)" = "scott", "Choose my own" = "fixed")),
          conditionalPanel("input.bw == 'fixed'",
                           sliderInput("sigma", "Bandwidth σ (km)", 2, 100, 15, step = 1))
        ),
        conditionalPanel(
          "input.method == 'ripley'",
          sliderInput("rmax", "Maximum distance r (km)", 5, 150, 50, step = 5),
          selectInput("nsim", "CSR simulations",
                      c("19 (quick)" = 19, "39 (5% test)" = 39, "99 (slow)" = 99), selected = 39),
          sliderInput("maxpts", "Max. points used", 500, 5000, 2000, step = 500)
        ),

        actionButton("run", "Run analysis", icon = icon("fire"), class = "btn-fire w-100"),
        tags$div(
          class = "sidebar-downloads",
          downloadButton("dl_results", "Results", class = "btn-outline-secondary btn-sm"),
          downloadButton("dl_points", "Detections", class = "btn-outline-secondary btn-sm")
        ),
        tags$p(class = "caveat",
               "Each point is a satellite-detected active-fire pixel (VIIRS 375 m), ",
               "not a verified wildfire or burned area.")
      ),

      layout_columns(
        fill = FALSE,
        value_box("Fire detections", textOutput("vb_n"),
                  showcase = bs_icon("fire"), tags$small(textOutput("vb_n_sub", inline = TRUE))),
        value_box("Nearest-neighbour ratio R", textOutput("vb_r"),
                  showcase = bs_icon("bullseye"), tags$small(textOutput("vb_r_sub", inline = TRUE))),
        value_box("Pattern classification", textOutput("vb_pattern"),
                  showcase = bs_icon("diagram-3"), tags$small(textOutput("vb_pattern_sub", inline = TRUE)))
      ),

      layout_columns(
        col_widths = c(7, 5),
        card(full_screen = TRUE, min_height = 460,
             card_header(textOutput("map_title", inline = TRUE)),
             leafletOutput("map", height = 500)),
        card(full_screen = TRUE, min_height = 460,
             card_header("Results & interpretation"),
             uiOutput("results"))
      )
    )
  ),

  # ---------------- Compare seasons ----------------
  nav_panel(
    "Compare seasons", icon = bs_icon("calendar3"),
    layout_sidebar(
      sidebar = sidebar(
        width = 300, title = "Season comparison",
        selectInput("s_region", "Region", REGIONS),
        selectInput("s_year", "Year", year_choices, selected = "all"),
        actionButton("s_run", "Compare seasons", icon = icon("fire"), class = "btn-fire w-100"),
        tags$p(class = "caveat",
               "Southern-hemisphere seasons: Summer Dec-Feb, Autumn Mar-May, ",
               "Winter Jun-Aug, Spring Sep-Nov. Each season is tested with the ",
               "nearest-neighbour method against the same study area.")
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(card_header("Detections per month"), plotOutput("monthly_plot", height = "420px")),
        card(card_header("Pattern by season"),
             uiOutput("season_text"), tableOutput("season_table"))
      )
    )
  ),

  # ---------------- About ----------------
  nav_panel(
    "About & methods", icon = bs_icon("info-circle"),
    layout_columns(
      col_widths = c(8, 4),
      card(card_body(class = "about", markdown(about_md()))),
      card(card_header("Data cleaning log"), tableOutput("clean_log"),
           tags$p(class = "caveat", textOutput("data_stamp", inline = TRUE)))
    )
  ),

  nav_spacer(),
  nav_item(tags$span(class = "nav-note", "NASA FIRMS · VIIRS S-NPP 375 m"))
)

# ---- server -----------------------------------------------------------------
server <- function(input, output, session) {

  filtered    <- reactive(filter_fires(input$region, input$year, input$month))
  current_key <- reactive(paste(input$region, input$year, input$month))
  res         <- reactiveVal(NULL)

  # ---------- value boxes ----------
  output$vb_n     <- renderText(fmt_num(nrow(filtered())))
  output$vb_n_sub <- renderText(paste(input$region, "·", period_label(input$year, input$month)))

  is_stale <- reactive(!is.null(res()) && res()$key != current_key())
  output$vb_r <- renderText(if (is.null(res())) "-" else fmt_num(res()$nn$R, 2))
  output$vb_r_sub <- renderText({
    if (is.null(res())) "Run an analysis to calculate"
    else if (is_stale()) "Settings changed - run again"
    else "< 1 clustered · 1 random · > 1 dispersed"
  })
  output$vb_pattern <- renderText(if (is.null(res())) "-" else res()$verdict)
  output$vb_pattern_sub <- renderText({
    if (is.null(res())) "Calculated live"
    else if (is_stale()) "Settings changed - run again"
    else paste("from", res()$verdict_source)
  })

  output$map_title <- renderText(
    paste(input$region, "·", period_label(input$year, input$month), "· fire detections"))

  # ---------- map ----------
  output$map <- renderLeaflet({
    leaflet(options = leafletOptions(preferCanvas = TRUE)) |>
      addProviderTiles(providers$CartoDB.Positron,  group = "Light") |>
      addProviderTiles(providers$CartoDB.DarkMatter, group = "Dark") |>
      addPolygons(data = provinces, color = FIRE$smoke, weight = 1, fillOpacity = 0,
                  label = ~province, group = "Provinces") |>
      addLayersControl(baseGroups = c("Light", "Dark"),
                       overlayGroups = c("Detections", "Density", "Provinces"),
                       options = layersControlOptions(collapsed = TRUE)) |>
      fitBounds(16.4, -34.9, 32.9, -22.1)
  })

  map_ready <- reactiveVal(FALSE)
  observeEvent(input$map_zoom, map_ready(TRUE), once = TRUE)

  # Zoom to and outline the selected region
  observe({
    req(map_ready())
    p <- leafletProxy("map") |> clearGroup("Selected")
    if (input$region == "South Africa") {
      p |> fitBounds(16.4, -34.9, 32.9, -22.1)
    } else {
      sel <- provinces[provinces$province == input$region, ]
      bb  <- as.numeric(st_bbox(sel))
      p |> addPolylines(data = sel, color = FIRE$ember, weight = 2.5, group = "Selected") |>
        fitBounds(bb[1], bb[2], bb[3], bb[4])
    }
  })

  # Detections as small ember dots (a random sample if there are very many)
  observe({
    req(map_ready())
    d <- filtered()
    if (nrow(d) > MAX_MAP_POINTS) {
      set.seed(1)
      d <- d[sample.int(nrow(d), MAX_MAP_POINTS), ]
    }
    p <- leafletProxy("map") |> clearGroup("Detections") |>
      clearGroup("Density") |> removeControl("kde_legend")      # old density is now stale
    if (nrow(d))
      p |> addCircleMarkers(data = d, lng = ~lon, lat = ~lat, radius = 2.5, stroke = FALSE,
                            fillColor = FIRE$ember, fillOpacity = 0.55, group = "Detections",
                            label = ~paste0(format(date, "%d %b %Y"), " · FRP ", frp, " MW"))
  })

  # ---------- run analysis ----------
  observeEvent(input$run, {
    d <- filtered()
    if (nrow(d) < 10) {
      showNotification("Fewer than 10 detections in this selection - choose a wider period or region.",
                       type = "warning")
      return()
    }
    region <- input$region
    label  <- paste(region, "(", period_label(input$year, input$month), ")")
    X <- make_ppp(d, windows[[region]])

    out <- tryCatch(withProgress(message = "Analysing the fire pattern...", value = 0.1, {
      o <- list(method = input$method, key = current_key(), label = label,
                region = region, period = period_label(input$year, input$month))
      o$nn <- run_nn(X)
      o$verdict <- o$nn$pattern
      o$verdict_source <- "nearest-neighbour test"
      incProgress(0.3)

      if (o$method == "kde") {
        o$kde <- run_kde(X, bw_method = input$bw, sigma_km = input$sigma)
        o$kde$peak_province <- province_at(o$kde$peak_lon, o$kde$peak_lat, provinces)
        o$text  <- interpret_kde(o$kde, label)
        o$table <- kde_table(o$kde)
      } else if (o$method == "nn") {
        o$text  <- interpret_nn(o$nn, label)
        o$table <- nn_table(o$nn)
      } else {
        incProgress(0.1, detail = "Simulating random patterns...")
        o$ripley <- run_ripley(X, rmax = input$rmax, nsim = as.integer(input$nsim),
                               max_points = input$maxpts)
        o$verdict <- o$ripley$verdict
        o$verdict_source <- "Ripley's K (global test)"
        o$text  <- interpret_ripley(o$ripley, label)
        o$table <- ripley_table(o$ripley)
      }
      incProgress(0.5)
      o
    }), error = function(e) {
      showNotification(paste("Analysis failed:", conditionMessage(e)), type = "error", duration = 10)
      NULL
    })
    req(out)
    res(out)

    if (out$method == "kde") {
      r   <- out$kde$raster
      rng <- range(terra::values(r), na.rm = TRUE)
      pal <- colorNumeric(KDE_PALETTE, domain = c(rng[1], rng[2] * 1.05), na.color = "transparent")
      leafletProxy("map") |>
        clearGroup("Density") |> removeControl("kde_legend") |>
        addRasterImage(r, colors = pal, opacity = 0.75, group = "Density", project = TRUE) |>
        addLegend(pal = pal, values = rng, title = "Detections<br>per 100 km²",
                  position = "bottomright", layerId = "kde_legend", group = "Density")
    }
  })

  # ---------- results panel ----------
  output$results <- renderUI({
    r <- res()
    if (is.null(r)) {
      return(tags$div(class = "empty-state", flame_svg(44),
                      tags$h5("Choose your settings, then click Run analysis"),
                      tags$p("KDE shows where detections concentrate, nearest neighbour tests ",
                             "whether they cluster, and Ripley's K shows at what distances.")))
    }
    method_name <- names(method_choices)[method_choices == r$method]
    tagList(
      tags$div(class = "result-head",
               tags$span(class = "method-chip", method_name),
               tags$span(class = paste("pattern-badge", tolower(r$verdict)), r$verdict)),
      if (is_stale()) tags$div(class = "stale-note", "Settings have changed since this result - click Run analysis to update."),
      tags$div(class = "callout", r$text),
      if (r$method == "ripley") tagList(
        radioButtons("ktype", NULL, c("L(r) - r (easier to read)" = "L", "K(r)" = "K"), inline = TRUE),
        plotOutput("ripley_plot", height = "300px")
      ),
      if (r$method == "nn") plotOutput("nn_plot", height = "220px"),
      tableOutput("res_table"),
      if (r$method == "kde") tags$p(class = "caveat", "The density layer has been added to the map (toggle it under the layers icon).")
    )
  })

  output$res_table   <- renderTable(res()$table, striped = TRUE, hover = TRUE, width = "100%")
  output$nn_plot     <- renderPlot(plot_nn_hist(res()$nn), res = 96)
  output$ripley_plot <- renderPlot({
    req(res()$ripley)
    plot_ripley(res()$ripley, type = if (is.null(input$ktype)) "L" else input$ktype)
  }, res = 96)

  # ---------- downloads ----------
  slug <- reactive(gsub("[^A-Za-z0-9]+", "_", tolower(current_key())))

  output$dl_results <- downloadHandler(
    filename = function() paste0("firewatch_", if (is.null(res())) "none" else res()$method, "_", slug(), ".csv"),
    content  = function(file) {
      r <- res()
      out <- if (is.null(r)) {
        data.frame(note = "No analysis has been run yet.")
      } else if (r$method == "ripley") {
        cv <- r$ripley$curve
        data.frame(region = r$region, period = r$period, r_km = cv$r,
                   L_obs = cv$obs, L_csr = cv$theo, L_lo = cv$lo, L_hi = cv$hi,
                   K_obs = pi * cv$obs^2, verdict = r$verdict, dclf_p = r$ripley$p_global)
      } else {
        cbind(region = r$region, period = r$period, method = r$method, r$table)
      }
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$dl_points <- downloadHandler(
    filename = function() paste0("firewatch_detections_", slug(), ".csv"),
    content  = function(file) {
      d <- filtered()
      write.csv(d[, c("date", "lon", "lat", "province", "season", "frp", "daynight")],
                file, row.names = FALSE)
    }
  )

  # ---------- compare seasons ----------
  season_data <- eventReactive(input$s_run, {
    d <- filter_fires(input$s_region, input$s_year, "all")
    list(d = d, region = input$s_region, year = input$s_year,
         table = compare_seasons(d, windows[[input$s_region]]))
  }, ignoreNULL = FALSE)

  output$monthly_plot <- renderPlot({
    s <- season_data()
    validate(need(nrow(s$d) > 0, "No detections for this selection."))
    plot_monthly(s$d)
  }, res = 96)

  output$season_table <- renderTable(season_data()$table, striped = TRUE, width = "100%")

  output$season_text <- renderUI({
    s <- season_data(); tb <- s$table
    req(sum(tb$Detections) > 0)
    top <- tb[which.max(tb$Detections), ]
    tags$div(class = "callout", sprintf(
      "In %s (%s), %s accounts for %s%% of all fire detections. %s",
      s$region, if (s$year == "all") "2024-2025" else s$year, top$Season, top$`Share (%)`,
      "Compare the NN ratios: values further below 1 mean detections in that season are more tightly clustered."))
  })

  # ---------- about ----------
  output$clean_log  <- renderTable(setNames(meta$cleaning, c("Step", "Rows")), digits = 0, width = "100%")
  output$data_stamp <- renderText(sprintf("Source: %s · %s to %s · prepared %s",
                                          paste(meta$source, collapse = ", "),
                                          meta$date_range[1], meta$date_range[2],
                                          format(meta$created, "%d %b %Y")))
}

shinyApp(ui, server)
