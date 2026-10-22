# Text for the "About & methods" tab (kept here so app.R stays readable)

about_md <- function() {
"
### FIREWATCH SA
**An Interactive Spatial Statistical Analysis of Satellite-Detected Fire Patterns in South Africa (2024-2025)**

The map shows where satellites detected active fires. The statistics test whether
the pattern we *see* is statistically meaningful, using three point-pattern methods.

#### Data
* **Fire detections:** NASA FIRMS, VIIRS S-NPP 375 m standard product (Collection 2,
  file `fire_archive_SV-C2_*`). One sensor is used throughout for consistency.
* **Cleaning:** low-confidence detections removed; only `type = 0` (presumed vegetation
  fire) kept, which removes industrial heat sources; exact duplicates removed.
* **Boundaries:** GADM v4.1 provinces. Analysis is done in South African Albers Equal
  Area (km), so densities are per km² and distances are in km.

#### Methods
| Method | Question | How to read it |
|---|---|---|
| **Kernel density (KDE)** | Where are detections concentrated? | A smooth surface of detections per 100 km². The bandwidth σ sets how much it is smoothed. |
| **Nearest neighbour (Clark-Evans)** | Clustered, random or dispersed? | R = observed ÷ expected mean nearest-neighbour distance. R < 1 clustered, R ≈ 1 random, R > 1 dispersed; tested with a z-test. |
| **Ripley's K / L function** | At what distances does clustering occur? | L(r) - r above the grey envelope of random (CSR) simulations means clustering at distance r; below means regular spacing. A global DCLF test gives one overall p-value. |

#### Limitations
* A detection is a 375 m pixel containing an active fire at overpass time - **not** a verified
  wildfire, and not burned area. One large fire can produce many detections.
* The tests compare against **complete spatial randomness** with constant intensity. Fire
  is naturally more likely in some landscapes (e.g. grasslands and savanna), so some
  'clustering' reflects land cover and not fire behaviour alone.
* The Ripley envelope is **pointwise**: crossing it at one distance is not by itself a
  formal test - use the global DCLF p-value.
* With very large samples, nearest-neighbour p-values are almost always tiny; the size of
  R is more informative than its p-value.
* Edge effects are handled with Ripley's border correction and KDE edge correction.

#### Acknowledgement
We acknowledge the use of data and imagery from NASA's Fire Information for Resource
Management System (FIRMS), part of NASA's Earth Science Data and Information System
(ESDIS). Boundaries: GADM.
"
}
