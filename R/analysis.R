# =============================================================================
# FIREWATCH SA - spatial point-pattern methods
#   1. Kernel Density Estimation (KDE)
#   2. Nearest-neighbour analysis (Clark-Evans ratio + z-test)
#   3. Ripley's K / L function with a CSR simulation envelope + global DCLF test
# All coordinates are in km (South African Albers), so areas are km² and
# distances are km.
# =============================================================================

# Build a spatstat point pattern from the cleaned detections table.
# Points were already checked to lie inside their province window during
# data preparation, so check = FALSE keeps this fast.
make_ppp <- function(df, W) {
  spatstat.geom::ppp(df$x, df$y, window = W, check = FALSE)
}

# -----------------------------------------------------------------------------
# 1. KERNEL DENSITY ESTIMATION
# -----------------------------------------------------------------------------
run_kde <- function(X, bw_method = c("scott", "fixed"), sigma_km = 15, dimyx = 256) {
  bw_method <- match.arg(bw_method)

  sigma <- if (bw_method == "scott") {
    s <- spatstat.explore::bw.scott(X)
    if (length(s) == 2) sqrt(prod(s)) else as.numeric(s)   # one isotropic value
  } else {
    sigma_km
  }

  # Gaussian kernel, edge-corrected, on a pixel grid. Units: detections per km²
  D <- spatstat.explore::density.ppp(X, sigma = sigma, dimyx = dimyx,
                                     edge = TRUE, positive = TRUE)

  v <- D$v * 100                      # -> detections per 100 km²
  vals <- v[!is.na(v)]

  # Concentration: share of total intensity found in the densest 10% of the area
  sorted   <- sort(vals, decreasing = TRUE)
  top_n    <- ceiling(0.10 * length(sorted))
  top10    <- sum(sorted[seq_len(top_n)]) / sum(sorted)

  # Peak location (pixel with the highest intensity)
  ij   <- which(v == max(vals), arr.ind = TRUE)[1, ]
  peak <- sf::st_sfc(sf::st_point(c(D$xcol[ij[2]], D$yrow[ij[1]])), crs = ALBERS_SA) |>
    sf::st_transform(4326) |>
    sf::st_coordinates()

  # Convert to a terra raster for Leaflet (row 1 must be the NORTHERN edge)
  r <- terra::rast(nrows = nrow(v), ncols = ncol(v),
                   xmin = D$xrange[1], xmax = D$xrange[2],
                   ymin = D$yrange[1], ymax = D$yrange[2],
                   crs = ALBERS_SA)
  terra::values(r) <- as.vector(t(v[nrow(v):1, ]))
  # Hide near-zero background so the map isn't tinted everywhere
  r[r < 0.02 * max(vals)] <- NA

  list(
    n            = spatstat.geom::npoints(X),
    sigma        = sigma,
    bw_method    = bw_method,
    raster       = r,
    peak_lon     = peak[1, "X"],
    peak_lat     = peak[1, "Y"],
    peak_density = max(vals),
    mean_density = mean(vals),
    top10_share  = top10
  )
}

interpret_kde <- function(k, region) {
  level <- if (k$top10_share >= 0.5) "highly concentrated" else
           if (k$top10_share >= 0.25) "moderately concentrated" else "fairly evenly spread"
  sprintf(paste(
    "Fire detections in %s are %s: the densest 10%% of the area holds %s%% of the",
    "estimated fire intensity (an even spread would give 10%%). The hotspot peaks near",
    "%.2f°E, %.2f°S%s at about %s detections per 100 km², compared with an average of %s.",
    "The bandwidth (σ = %s km) controls the smoothing: smaller values reveal local hotspots,",
    "larger values show the regional trend."),
    region, level, fmt_num(100 * k$top10_share, 0),
    k$peak_lon, abs(k$peak_lat),
    if (!is.null(k$peak_province) && !is.na(k$peak_province)) paste0(" (", k$peak_province, ")") else "",
    fmt_num(k$peak_density, 1), fmt_num(k$mean_density, 1), fmt_num(k$sigma, 1))
}

kde_table <- function(k) {
  data.frame(
    Statistic = c("Detections", "Bandwidth σ", "Peak intensity", "Mean intensity",
                  "Intensity in densest 10% of area", "Peak location"),
    Value = c(fmt_num(k$n),
              paste(fmt_num(k$sigma, 1), "km", if (k$bw_method == "scott") "(Scott's rule)" else "(fixed)"),
              paste(fmt_num(k$peak_density, 1), "per 100 km²"),
              paste(fmt_num(k$mean_density, 1), "per 100 km²"),
              paste0(fmt_num(100 * k$top10_share, 1), "%"),
              sprintf("%.2f°E, %.2f°S", k$peak_lon, abs(k$peak_lat))),
    check.names = FALSE
  )
}

# -----------------------------------------------------------------------------
# 2. NEAREST-NEIGHBOUR ANALYSIS (Clark & Evans, 1954)
#    R = observed mean NN distance / expected mean NN distance under CSR
#    expected = 0.5 / sqrt(lambda),  SE = 0.26136 / sqrt(n * lambda)
# -----------------------------------------------------------------------------
run_nn <- function(X, alpha = 0.05) {
  n <- spatstat.geom::npoints(X)
  if (n < 3) return(NULL)
  A      <- spatstat.geom::area.owin(spatstat.geom::Window(X))
  lambda <- n / A
  d      <- spatstat.geom::nndist(X)

  observed <- mean(d)
  expected <- 0.5 / sqrt(lambda)
  se       <- 0.26136 / sqrt(n * lambda)
  R        <- observed / expected
  z        <- (observed - expected) / se
  p        <- 2 * stats::pnorm(-abs(z))
  pattern  <- if (p >= alpha) "Random" else if (R < 1) "Clustered" else "Dispersed"

  list(n = n, area_km2 = A, intensity = lambda, observed = observed,
       expected = expected, R = R, z = z, p = p, pattern = pattern, distances = d)
}

interpret_nn <- function(nn, region) {
  direction <- switch(nn$pattern,
    Clustered = "significantly CLUSTERED - detections sit closer together than chance would produce",
    Dispersed = "significantly DISPERSED - detections are spaced out more evenly than chance",
    Random    = "not significantly different from RANDOM (complete spatial randomness)")
  big_n <- if (nn$n > 10000)
    " Note: with this many detections even tiny departures from randomness are 'significant', so focus on how far R is from 1." else ""
  sprintf(paste(
    "R = %s. On average a detection's nearest neighbour is %s km away, versus %s km expected",
    "if the same %s detections were scattered at random over %s. The z-test (z = %s, p %s)",
    "shows the pattern is %s.%s"),
    fmt_num(nn$R, 2), fmt_num(nn$observed, 2), fmt_num(nn$expected, 2),
    fmt_num(nn$n), region, fmt_num(nn$z, 2),
    if (nn$p < 0.001) "< 0.001" else paste("=", fmt_p(nn$p)),
    direction, big_n)
}

nn_table <- function(nn) {
  data.frame(
    Statistic = c("Detections (n)", "Study area", "Intensity (λ)",
                  "Observed mean NN distance", "Expected under CSR",
                  "Clark-Evans ratio R", "z-score", "p-value", "Pattern"),
    Value = c(fmt_num(nn$n), paste(fmt_num(nn$area_km2), "km²"),
              paste(fmt_num(100 * nn$intensity, 2), "per 100 km²"),
              paste(fmt_num(nn$observed, 3), "km"), paste(fmt_num(nn$expected, 3), "km"),
              fmt_num(nn$R, 3), fmt_num(nn$z, 2), fmt_p(nn$p), nn$pattern),
    check.names = FALSE
  )
}

plot_nn_hist <- function(nn) {
  d <- nn$distances
  d <- d[d <= stats::quantile(d, 0.99)]          # trim the long tail for readability
  ggplot2::ggplot(data.frame(d = d), ggplot2::aes(d)) +
    ggplot2::geom_histogram(bins = 40, fill = FIRE$flame, colour = "white", linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = nn$expected, linetype = "dashed", colour = FIRE$smoke) +
    ggplot2::annotate("text", x = nn$expected, y = Inf, label = " expected under CSR",
                      hjust = 0, vjust = 1.5, size = 3.5, colour = FIRE$smoke) +
    ggplot2::labs(x = "Distance to nearest detection (km)", y = "Detections") +
    ggplot2::theme_minimal(base_size = 12)
}

# -----------------------------------------------------------------------------
# 3. RIPLEY'S K / L FUNCTION
#    L(r) = sqrt(K(r)/pi);  under CSR L(r) = r, so L(r) - r = 0.
#    Above the envelope -> clustering at that distance; below -> regularity.
#    Random thinning preserves K, so a random subsample gives the same K(r)
#    as the full pattern (only with more noise) - this keeps it fast.
# -----------------------------------------------------------------------------
run_ripley <- function(X, rmax = 50, nsim = 39, max_points = 2000, seed = 2026) {
  n_all <- spatstat.geom::npoints(X)
  if (n_all > max_points) {
    set.seed(seed)
    X <- X[sample.int(n_all, max_points)]
  }
  r <- seq(0, rmax, length.out = 101)

  env <- spatstat.explore::envelope(
    X, fun = spatstat.explore::Lest, nsim = nsim, r = r,
    correction = "border", savefuns = TRUE, verbose = FALSE
  )
  curve <- as.data.frame(env)[, c("r", "obs", "theo", "lo", "hi")]

  # Global DCLF test re-uses the same simulations (one p-value for all r)
  p_global <- tryCatch(
    spatstat.explore::dclf.test(env, rinterval = c(0, rmax))$p.value,
    error = function(e) NA_real_)

  ok    <- curve$r > 0 & !is.na(curve$obs)
  above <- ok & curve$obs > curve$hi
  below <- ok & curve$obs < curve$lo
  # Direction of the departure: how far the curve goes above vs below the envelope
  excess_up   <- sum(pmax(curve$obs - curve$hi, 0)[ok])
  excess_down <- sum(pmax(curve$lo - curve$obs, 0)[ok])
  dev <- if (excess_up + excess_down > 0) excess_up - excess_down else mean((curve$obs - curve$theo)[ok])

  verdict <- if (!is.na(p_global)) {
    if (p_global >= 0.05) "Random" else if (dev > 0) "Clustered" else "Dispersed"
  } else {                                   # fallback if the global test fails
    if (mean(above[ok]) > 0.2) "Clustered" else if (mean(below[ok]) > 0.2) "Dispersed" else "Random"
  }

  diff_lr <- curve$obs - curve$r
  peak_r  <- if (any(above)) curve$r[above][which.max(diff_lr[above])] else NA

  list(curve = curve, n_used = spatstat.geom::npoints(X), n_all = n_all,
       nsim = nsim, rmax = rmax, p_global = p_global, verdict = verdict,
       clustered_at = describe_runs(curve$r, above),
       regular_at   = describe_runs(curve$r, below),
       peak_r = peak_r)
}

# Turn a TRUE/FALSE vector along r into readable distance ranges, e.g. "2-35 km"
describe_runs <- function(r, flag) {
  flag[is.na(flag)] <- FALSE
  if (!any(flag)) return(character(0))
  rl     <- rle(flag)
  ends   <- cumsum(rl$lengths)
  starts <- ends - rl$lengths + 1
  keep   <- rl$values
  sprintf("%s-%s km", fmt_num(r[starts[keep]], 1), fmt_num(r[ends[keep]], 1))
}

interpret_ripley <- function(k, region) {
  sample_note <- if (k$n_used < k$n_all)
    sprintf(" A random sample of %s of the %s detections was used (random thinning does not change K(r)).",
            fmt_num(k$n_used), fmt_num(k$n_all)) else ""
  where <- if (length(k$clustered_at))
    sprintf("The observed curve rises above the CSR envelope at %s, so detections are clustered at those distances%s.",
            paste(k$clustered_at, collapse = ", "),
            if (!is.na(k$peak_r)) sprintf(", most strongly around r ≈ %s km", fmt_num(k$peak_r, 0)) else "")
  else "The observed curve does not rise above the CSR envelope at any distance tested."
  reg <- if (length(k$regular_at))
    sprintf(" It falls below the envelope (regular spacing) at %s.", paste(k$regular_at, collapse = ", ")) else ""
  sprintf("%s%s Overall (global DCLF test, %d simulations, p %s) the pattern in %s is classed as %s.%s",
          where, reg, k$nsim,
          if (is.na(k$p_global)) "not available" else if (k$p_global < 0.001) "< 0.001" else paste("=", fmt_p(k$p_global)),
          region, toupper(k$verdict), sample_note)
}

ripley_table <- function(k) {
  data.frame(
    Statistic = c("Points used", "CSR simulations", "Distances tested",
                  "Clustering at", "Regularity at", "Strongest clustering",
                  "Global DCLF p-value", "Verdict"),
    Value = c(paste(fmt_num(k$n_used), "of", fmt_num(k$n_all)), k$nsim,
              paste0("0-", k$rmax, " km"),
              if (length(k$clustered_at)) paste(k$clustered_at, collapse = ", ") else "none",
              if (length(k$regular_at))   paste(k$regular_at,   collapse = ", ") else "none",
              if (is.na(k$peak_r)) "-" else paste("r ≈", fmt_num(k$peak_r, 0), "km"),
              fmt_p(k$p_global), k$verdict),
    check.names = FALSE
  )
}

plot_ripley <- function(k, type = c("L", "K")) {
  type <- match.arg(type)
  cv <- k$curve
  if (type == "L") {
    df <- data.frame(r = cv$r, obs = cv$obs - cv$r, theo = cv$theo - cv$r,
                     lo = cv$lo - cv$r, hi = cv$hi - cv$r)
    ylab <- "L(r) - r"
  } else {                                  # K = pi * L^2 (monotone, so envelope transfers)
    df <- data.frame(r = cv$r, obs = pi * cv$obs^2, theo = pi * cv$theo^2,
                     lo = pi * cv$lo^2, hi = pi * cv$hi^2)
    ylab <- "K(r)  (km²)"
  }
  ggplot2::ggplot(df, ggplot2::aes(r)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = "CSR envelope"), alpha = 0.6) +
    ggplot2::geom_line(ggplot2::aes(y = theo, colour = "Expected under CSR"), linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(y = obs, colour = "Observed"), linewidth = 1.1) +
    ggplot2::scale_fill_manual(NULL, values = c("CSR envelope" = "#E9E2DA")) +
    ggplot2::scale_colour_manual(NULL, values = c("Observed" = FIRE$ember, "Expected under CSR" = FIRE$smoke)) +
    ggplot2::labs(x = "Distance r (km)", y = ylab) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

# -----------------------------------------------------------------------------
# SEASON COMPARISON (nearest-neighbour per season + monthly counts)
# -----------------------------------------------------------------------------
compare_seasons <- function(df, W) {
  total <- nrow(df)
  rows <- lapply(SEASON_LEVELS, function(s) {
    sub <- df[df$season == s, ]
    nn  <- if (nrow(sub) >= 3) run_nn(make_ppp(sub, W)) else NULL
    data.frame(
      Season       = s,
      Detections   = nrow(sub),
      `Share (%)`  = if (total) round(100 * nrow(sub) / total, 1) else NA,
      `NN ratio R` = if (is.null(nn)) NA else round(nn$R, 3),
      `z-score`    = if (is.null(nn)) NA else round(nn$z, 2),
      `p-value`    = if (is.null(nn)) "-" else fmt_p(nn$p),
      Pattern      = if (is.null(nn)) "Too few points" else nn$pattern,
      check.names  = FALSE
    )
  })
  do.call(rbind, rows)
}

plot_monthly <- function(df) {
  counts <- as.data.frame(table(year = df$year, month = df$month))
  counts$month  <- factor(month.abb[as.integer(as.character(counts$month))], levels = month.abb)
  counts$season <- month_to_season(match(counts$month, month.abb))
  ggplot2::ggplot(counts, ggplot2::aes(month, Freq, fill = season)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_wrap(~ year, ncol = 1) +
    ggplot2::scale_fill_manual(NULL, values = SEASON_COLOURS, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = function(x) format(x, big.mark = " ", scientific = FALSE)) +
    ggplot2::labs(x = NULL, y = "Fire detections") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.major.x = ggplot2::element_blank())
}

# Which province contains a lon/lat point?
province_at <- function(lon, lat, provinces) {
  old <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old)))
  pt <- sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326)
  i  <- tryCatch(sf::st_intersects(pt, provinces)[[1]], error = function(e) integer(0))
  if (length(i)) as.character(provinces$province[i[1]]) else NA_character_
}
