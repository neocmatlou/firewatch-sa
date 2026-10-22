# =============================================================================
# FIREWATCH SA - method validation (member 2)
# Simulates three patterns whose answer we KNOW, then checks that the
# nearest-neighbour test and Ripley's K classify each one correctly.
# Run from the project root:  source("tests/validate_methods.R")
# =============================================================================

source("R/constants.R")
source("R/analysis.R")
suppressPackageStartupMessages({
  library(spatstat.geom); library(spatstat.random); library(spatstat.explore)
})

set.seed(2026)
W <- owin(c(0, 100), c(0, 100))                                        # 100 km x 100 km

patterns <- list(
  Random    = rpoispp(0.05, win = W),                                  # CSR, ~500 points
  Clustered = rThomas(kappa = 0.005, scale = 3, mu = 10, win = W),     # ~50 clusters of ~10
  Dispersed = rSSI(r = 3, n = 400, win = W)                            # inhibition: no pair < 3 km
)

results <- do.call(rbind, lapply(names(patterns), function(truth) {
  X  <- patterns[[truth]]
  nn <- run_nn(X)
  rk <- run_ripley(X, rmax = 20, nsim = 39)
  data.frame(
    Truth          = truth,
    n              = npoints(X),
    `NN R`         = round(nn$R, 3),
    `NN p`         = fmt_p(nn$p),
    `NN verdict`   = nn$pattern,
    `Ripley p`     = fmt_p(rk$p_global),
    `Ripley verdict` = rk$verdict,
    Pass           = ifelse(nn$pattern == truth & rk$verdict == truth, "PASS", "FAIL"),
    check.names = FALSE
  )
}))

print(results, row.names = FALSE)

# Picture for the presentation
op <- par(mfrow = c(1, 3), mar = c(1, 1, 3, 1))
for (nm in names(patterns)) plot(patterns[[nm]], main = nm, pch = 16, cex = 0.4, cols = FIRE$ember)
par(op)

if (all(results$Pass == "PASS")) {
  message("\nAll three patterns classified correctly by both methods.")
} else {
  message("\nAt least one check failed - with random data this can happen ~5% of the time ",
          "for the CSR pattern. Re-run with another seed; if it fails again, check R/analysis.R.")
}
