# =============================================================================
# FIREWATCH SA - data preparation (member 1)
# Input : data-raw/fire_archive_SV-C2_828855.csv  (or .xlsx)  from NASA FIRMS
#         South African province boundaries (GADM, downloaded automatically)
# Output: data/fires.rds      cleaned detections (lon/lat + Albers km x/y)
#         data/provinces.rds  simplified province polygons (WGS84, for the map)
#         data/windows.rds    spatstat observation windows (9 provinces + SA)
#         data/meta.rds       years, cleaning log, source info
#
# Run from the project root:  source("scripts/01_prepare_data.R")
# =============================================================================

source("R/constants.R")
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(readxl); library(lubridate)
  library(sf); library(spatstat.geom); library(geodata)
})

# ---- settings ---------------------------------------------------------------
TEST_MODE   <- TRUE        # TRUE = quick check on Limpopo, August 2024 (nothing saved)
YEARS       <- 2024:2025
SIMPLIFY_KM <- 1           # boundary simplification tolerance (km)

dir.create("data", showWarnings = FALSE)
dir.create("data-raw", showWarnings = FALSE)
log_step <- function(log, step, n) { message(sprintf("%-45s %10s rows", step, format(n, big.mark = " "))); rbind(log, data.frame(step = step, rows = n)) }
clean_log <- data.frame(step = character(), rows = integer())

# ---- 1. find and read the FIRMS file(s) ------------------------------------
raw_files <- list.files("data-raw", pattern = "^fire_archive.*\\.(csv|xlsx|xls)$",
                        full.names = TRUE, ignore.case = TRUE)
if (!length(raw_files))
  stop("No FIRMS file found. Put fire_archive_SV-C2_828855.csv (or .xlsx) in data-raw/")
if (length(list.files("data-raw", pattern = "^fire_nrt", ignore.case = TRUE)))
  warning("A fire_nrt_* file is present - it is IGNORED. Do not mix near-real-time and standard (archive) data.")

read_firms <- function(path) {
  message("Reading ", basename(path))
  ext <- tolower(tools::file_ext(path))
  df <- if (ext == "csv") {
    readr::read_csv(path, col_types = readr::cols(.default = "c"), progress = FALSE)
  } else {
    readxl::read_excel(path, col_types = "text")
  }
  names(df) <- tolower(trimws(names(df)))
  df
}

# FIRMS dates are "2024-08-15", but Excel may have changed them to 15/08/2024
# or to a serial number like 45519. Handle all three.
parse_firms_date <- function(x) {
  x <- trimws(x)
  out <- as.Date(rep(NA_character_, length(x)))
  serial <- grepl("^[0-9]{5}(\\.0+)?$", x)
  out[serial]  <- as.Date(as.numeric(x[serial]), origin = "1899-12-30")
  out[!serial] <- as.Date(lubridate::parse_date_time(x[!serial], orders = c("Ymd", "dmY", "mdY", "Ymd HMS"), quiet = TRUE))
  out
}

raw <- bind_rows(lapply(raw_files, read_firms))
needed <- c("latitude", "longitude", "acq_date", "confidence")
if (!all(needed %in% names(raw)))
  stop("Missing columns: ", paste(setdiff(needed, names(raw)), collapse = ", "),
       "\nColumns found: ", paste(names(raw), collapse = ", "))
for (col in c("acq_time", "type", "frp", "daynight", "satellite", "instrument"))
  if (!col %in% names(raw)) raw[[col]] <- NA_character_

clean_log <- log_step(clean_log, "Raw detections read", nrow(raw))
if (nrow(raw) %in% c(1048575, 1048576))
  warning("Exactly Excel's row limit - the file was probably truncated by Excel. Use the original CSV from FIRMS.")

fires <- raw |>
  transmute(
    lat       = as.numeric(latitude),
    lon       = as.numeric(longitude),
    date      = parse_firms_date(acq_date),
    acq_time  = acq_time,
    conf_raw  = tolower(trimws(confidence)),
    type      = suppressWarnings(as.integer(type)),
    frp       = suppressWarnings(as.numeric(frp)),
    daynight  = toupper(daynight),
    satellite = satellite
  )
if (mean(is.na(fires$date)) > 0.01)
  stop("More than 1% of dates could not be read. Example values: ",
       paste(head(raw$acq_date[is.na(fires$date)], 3), collapse = ", "))
fires <- fires |> filter(!is.na(lat), !is.na(lon), !is.na(date))
message("Instrument / satellite in file: ",
        paste(unique(na.omit(raw$instrument)), collapse = ", "), " / ",
        paste(unique(na.omit(raw$satellite)),  collapse = ", "))

# ---- 2. clean ---------------------------------------------------------------
fires <- fires |> filter(year(date) %in% YEARS)
clean_log <- log_step(clean_log, paste("Within", min(YEARS), "-", max(YEARS)), nrow(fires))

if (TEST_MODE) {
  fires <- fires |> filter(year(date) == 2024, month(date) == 8)
  clean_log <- log_step(clean_log, "TEST MODE: August 2024 only", nrow(fires))
}

# Confidence: VIIRS uses l / n / h (low, nominal, high); MODIS uses 0-100
conf_num <- suppressWarnings(as.numeric(fires$conf_raw))
keep_conf <- ifelse(!is.na(conf_num), conf_num >= 30, substr(fires$conf_raw, 1, 1) %in% c("n", "h"))
fires <- fires[keep_conf, ]
clean_log <- log_step(clean_log, "Removed low-confidence detections", nrow(fires))

# Type: 0 = presumed vegetation fire, 1 = volcano, 2 = other static land source
#       (e.g. industrial heat), 3 = offshore. Keep vegetation fires only.
if (all(is.na(fires$type))) {
  warning("No 'type' column (NRT data?). Industrial hot spots cannot be removed.")
} else {
  fires <- fires |> filter(type == 0)
  clean_log <- log_step(clean_log, "Kept vegetation fires only (type = 0)", nrow(fires))
}

fires <- fires |> distinct(lat, lon, date, acq_time, .keep_all = TRUE)
clean_log <- log_step(clean_log, "Removed exact duplicates", nrow(fires))

# ---- 3. province boundaries (GADM level 1) ---------------------------------
gpkg <- list.files("data-raw", pattern = "gadm41_ZAF.*\\.gpkg$", full.names = TRUE)
prov <- if (length(gpkg)) {
  st_read(gpkg[1], layer = "ADM_ADM_1", quiet = TRUE)
} else {
  message("Downloading GADM boundaries...")
  g <- geodata::gadm("ZAF", level = 1, path = "data-raw")
  if (is.null(g)) stop("GADM download failed. Download 'gadm41_ZAF.gpkg' manually from ",
                       "https://gadm.org/download_country.html (South Africa, GeoPackage) into data-raw/ and re-run.")
  st_as_sf(g)
}
prov <- prov |> select(province = NAME_1)
prov$province[prov$province == "Kwazulu-Natal"] <- "KwaZulu-Natal"   # spelling safety
prov$province[prov$province == "North-West"]    <- "North West"

# Drop the Prince Edward Islands (they would stretch the window into the ocean)
suppressMessages(sf_use_s2(FALSE))
prov <- suppressWarnings(st_crop(st_make_valid(prov), xmin = 16, ymin = -35.2, xmax = 33.2, ymax = -22))
suppressMessages(sf_use_s2(TRUE))

prov_km <- prov |>
  st_transform(ALBERS_SA) |>
  st_simplify(preserveTopology = TRUE, dTolerance = SIMPLIFY_KM) |>
  st_make_valid() |>
  (\(g) suppressWarnings(st_collection_extract(g, "POLYGON")))() |>
  group_by(province) |> summarise(.groups = "drop") |>
  arrange(province)

stopifnot(setequal(prov_km$province, PROVINCES))

# ---- 4. spatstat observation windows (km) ----------------------------------
to_owin <- function(g) {
  tryCatch(as.owin(g), error = function(e) as.owin(st_buffer(g, 0)))
}
windows <- lapply(seq_len(nrow(prov_km)), function(i) to_owin(st_geometry(prov_km)[i]))
names(windows) <- prov_km$province
windows[["South Africa"]] <- to_owin(st_union(prov_km))
message("Country window area: ", format(round(area.owin(windows[["South Africa"]])), big.mark = " "), " km²")

# ---- 5. project points, assign provinces -----------------------------------
pts <- st_as_sf(fires, coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  st_transform(ALBERS_SA)
pts <- st_join(pts, prov_km, join = st_intersects, left = FALSE)
xy  <- st_coordinates(pts)
fires <- st_drop_geometry(pts) |>
  mutate(x = xy[, 1], y = xy[, 2]) |>
  distinct(lat, lon, date, acq_time, .keep_all = TRUE)        # border points can match twice
clean_log <- log_step(clean_log, "Inside South Africa", nrow(fires))

# Keep only points that fall inside their (simplified) province window
inside <- logical(nrow(fires))
for (p in PROVINCES) {
  idx <- which(fires$province == p)
  if (length(idx)) inside[idx] <- inside.owin(fires$x[idx], fires$y[idx], windows[[p]])
}
fires <- fires[inside, ]
clean_log <- log_step(clean_log, "Inside simplified province windows", nrow(fires))

# ---- 6. final columns -------------------------------------------------------
fires <- fires |>
  transmute(
    lon = round(lon, 5), lat = round(lat, 5),
    x = round(x, 3), y = round(y, 3),                 # km, 1 m precision
    date, year = as.integer(year(date)), month = as.integer(month(date)),
    season = month_to_season(month),
    province = factor(province, levels = PROVINCES),
    frp = round(frp, 1),
    daynight = factor(daynight, levels = c("D", "N"), labels = c("Day", "Night"))
  ) |>
  arrange(date)

# ---- 7. sanity checks -------------------------------------------------------
message("\nDetections by province and year:")
print(table(fires$province, fires$year))
message("\nDetections by month:")
print(table(factor(month.abb[fires$month], levels = month.abb), fires$year))

if (TEST_MODE) {
  lim <- fires |> filter(province == "Limpopo")
  message("\nTEST MODE: Limpopo, August 2024 has ", nrow(lim), " detections.")
  if (nrow(lim) > 10) {
    source("R/analysis.R")
    nn <- run_nn(make_ppp(lim, windows[["Limpopo"]]))
    message("Clark-Evans R = ", round(nn$R, 3), "  (", nn$pattern, ")")
  }
  message("\nTest passed - set TEST_MODE <- FALSE and run again to build data/.")
} else {
  provinces <- st_transform(prov_km, 4326)
  meta <- list(
    years       = sort(unique(fires$year)),
    n_total     = nrow(fires),
    date_range  = range(fires$date),
    cleaning    = clean_log,
    source      = basename(raw_files),
    created     = Sys.time()
  )
  saveRDS(fires,     "data/fires.rds",     compress = "xz")
  saveRDS(provinces, "data/provinces.rds", compress = "xz")
  saveRDS(windows,   "data/windows.rds",   compress = "xz")
  saveRDS(meta,      "data/meta.rds")
  message("\nSaved data/*.rds  (", format(nrow(fires), big.mark = " "), " detections). ",
          "fires.rds is ", round(file.size("data/fires.rds") / 1e6, 1), " MB.")
  print(clean_log)
}
