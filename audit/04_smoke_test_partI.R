# Stage 1 smoke test: run the historical part I script VERBATIM (commit
# 74befc2) for the total canopy production results, patching only file paths.
# `ten plot trees.csv` is the Figshare copy (file 29128734), committed to
# data/ in f20319b.
# Run from the repository root:  Rscript audit/04_smoke_test_partI.R

treePath <- "data/ten plot trees.csv"

script <- readLines("scripts/historical/part I (updated with spp by plot).R")
run <- script[-(1:6)]  # drop setwd()
run <- sub('read.csv("ten plot trees.csv"', 'read.csv(treePath', run,
           fixed = TRUE)
run <- sub('read.csv("', 'read.csv("data/', run, fixed = TRUE)
run <- sub('read.csv("data/SSI Jan21.csv"',
           'read.csv("../ecophysio_traits/raw_data/SSI Jan21.csv"', run,
           fixed = TRUE)

pdf(NULL)
set.seed(1)
sink(tempfile())  # silence the many summary() prints
for (e in parse(text = run)) {
  res <- try(eval(e, envir = globalenv()), silent = TRUE)
  if (inherits(res, "try-error")) {
    sink(); cat("FAILED:", substr(deparse(e)[1], 1, 70), "\n")
    sink(tempfile())
  }
}
sink()

cat("\n== Index-based selections (fragile: depend on factor-level order) ==\n")
cat("lianas:", lianas, sep = "\n  ")
cat("spp.sel (n =", length(spp.sel), "):", spp.sel, sep = "\n  ")

cat("\n== Feb21 CSV reproduced? ==\n")
feb21 <- read.csv("data/litter.production w species-plot CAA Feb21.csv")
lc <- litter.consol[order(litter.consol$species, litter.consol$plot), ]
f <- feb21[order(feb21$species, feb21$plot), ]
cat("rows", nrow(lc), "vs", nrow(f), "\n")
for (v in c("litter.collection", "ten.ba", "duration", "litter.production",
            "SSI", "iv.wet", "CNratio"))
  cat(sprintf("  %-18s max abs diff %.3g\n", v,
              max(abs(lc[[v]] - f[[v]]), na.rm = TRUE)))

cat("\n== Total canopy production (paper: 768 +/- 48 g/m2/yr) ==\n")
cat("mean", mean(all), "SE", se(all), "\n")
cat("t-tests (paper: all t = 1.33 p = 0.221; no lianas t = 0.49 p = 0.639;",
    "per BA t = 6.81 p < 0.001)\n")
for (y in list(all, bolian, bolianba)) {
  tt <- t.test(y ~ pt, var.equal = TRUE)
  cat(sprintf("  t = %.2f df = %d p = %.4f\n", tt$statistic, tt$parameter,
              tt$p.value))
}
