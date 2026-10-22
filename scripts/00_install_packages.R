# Run once per computer:  source("scripts/00_install_packages.R")

pkgs <- c(
  # dashboard
  "shiny", "bslib", "bsicons", "leaflet", "ggplot2",
  # spatial
  "sf", "terra", "geodata", "spatstat.geom", "spatstat.explore", "spatstat.random",
  # data preparation
  "dplyr", "readr", "readxl", "lubridate",
  # deployment
  "rsconnect"
)

missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing) else message("All packages already installed.")

# Quick check that the spatial stack loads
invisible(lapply(c("sf", "terra", "spatstat.geom", "spatstat.explore", "leaflet"),
                 library, character.only = TRUE))
message("Ready. sf ", packageVersion("sf"), " | spatstat.explore ", packageVersion("spatstat.explore"),
        " | leaflet ", packageVersion("leaflet"))
