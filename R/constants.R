# =============================================================================
# FIREWATCH SA - shared constants and small helpers
# Shiny auto-loads every file in R/ when the app starts.
# The scripts in scripts/ and tests/ source this file explicitly.
# =============================================================================

# South African Albers Equal Area projection, in KILOMETRES.
# Equal-area + km units means densities are "per km²" and distances are "km".
ALBERS_SA <- paste(
  "+proj=aea +lat_1=-24 +lat_2=-33 +lat_0=0 +lon_0=25",
  "+x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
)

# Fire palette (used sparingly: accents, the density layer and charts)
FIRE <- list(
  ember    = "#D9480F",
  flame    = "#F08C00",
  amber    = "#FFB703",
  deep     = "#9D0208",
  charcoal = "#1F1A17",
  ash      = "#FBF8F4",
  smoke    = "#6B625C"
)

# Kernel-density colour ramp: pale ash -> amber -> ember -> deep red
KDE_PALETTE <- c("#FFF4C2", "#FFC857", "#F08C00", "#D9480F", "#9D0208", "#370617")

PROVINCES <- c(
  "Eastern Cape", "Free State", "Gauteng", "KwaZulu-Natal", "Limpopo",
  "Mpumalanga", "North West", "Northern Cape", "Western Cape"
)
REGIONS <- c("South Africa", PROVINCES)

# Southern-hemisphere meteorological seasons
SEASON_LEVELS <- c("Summer (Dec-Feb)", "Autumn (Mar-May)", "Winter (Jun-Aug)", "Spring (Sep-Nov)")
SEASON_COLOURS <- setNames(c("#D9480F", "#B08968", "#5C677D", "#F08C00"), SEASON_LEVELS)

month_to_season <- function(m) {
  s <- ifelse(m %in% c(12, 1, 2), SEASON_LEVELS[1],
       ifelse(m %in% 3:5,        SEASON_LEVELS[2],
       ifelse(m %in% 6:8,        SEASON_LEVELS[3], SEASON_LEVELS[4])))
  factor(s, levels = SEASON_LEVELS)
}

# Maximum number of detections drawn as dots on the map (statistics always use all)
MAX_MAP_POINTS <- 8000

# ---- formatting helpers -----------------------------------------------------
fmt_num <- function(x, digits = 0) {
  ifelse(is.na(x), "-", formatC(x, format = "f", digits = digits, big.mark = " "))
}

fmt_p <- function(p) {
  if (is.na(p)) return("-")
  if (p < 0.001) "< 0.001" else formatC(p, format = "f", digits = 3)
}

period_label <- function(year, month) {
  y <- if (year == "all") "2024-2025" else year
  m <- if (month == "all") "all months" else month.name[as.integer(month)]
  paste(m, y)
}

# Small inline flame logo (the one decorative "fire" element in the header)
flame_svg <- function(size = 22) {
  htmltools::HTML(sprintf(
    '<svg class="flame-logo" width="%d" height="%d" viewBox="0 0 24 24" aria-hidden="true">
       <defs><linearGradient id="fwFlame" x1="0" y1="1" x2="0" y2="0">
         <stop offset="0" stop-color="#9D0208"/><stop offset="0.55" stop-color="#D9480F"/>
         <stop offset="1" stop-color="#FFB703"/></linearGradient></defs>
       <path fill="url(#fwFlame)" d="M12 1.5c.9 3.1-.8 5.2-2.3 7.1C8.2 10.5 6.5 12.4 6.5 15.4a5.5 5.5 0 0 0 11 0c0-2.3-1.1-4-2.2-5.2-.1 1.7-.8 2.9-1.9 3.5.5-3.2-.4-6.6-1.4-8.2z"/>
     </svg>', size, size))
}
