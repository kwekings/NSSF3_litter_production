# Stage 1 smoke test: run the historical main script VERBATIM (commit 74befc2)
# up to the component models, patching only file paths, and compare the
# results that do not need piecewiseSEM with the published paper.
#
# Nothing here is the refactor; it only tests whether the historical code +
# historical inputs reproduce the published numbers.
# Run from the repository root:  Rscript audit/03_smoke_test_main.R

library(lme4)

script <- system(paste0('git show "74befc2:scripts/',
                        'leaf traits updated with spp by plot and iv wet.R"'),
                 intern = TRUE)

# Lines 4-226: data, PCA, merges, trait/SIV/canopy component models.
# Line 91 calls pairs.cor(), which is only defined at line 785, so it is
# skipped (a plot; no effect on results). Lines 1-2 are setwd().
run <- script[4:226]
run[91 - 3] <- "# pairs.cor() skipped"
run <- sub('read.csv("', 'read.csv("data/', run, fixed = TRUE)
# piecewiseSEM is not installed yet (version choice pending); nothing in
# lines 4-226 uses it
run <- sub("library(piecewiseSEM)", "# library(piecewiseSEM)", run,
           fixed = TRUE)

pdf(NULL)
# Evaluate expression by expression: models using SSI.stem fail because the
# Feb21 CSV has no SSI.stem column (see AGENTS.md). Log failures, carry on.
for (e in parse(text = run)) {
  res <- try(eval(e, envir = globalenv()), silent = TRUE)
  if (inherits(res, "try-error"))
    cat("FAILED:", substr(deparse(e)[1], 1, 70), "\n")
}

cat("\n== Sample sizes ==\n")
cat("species:", length(unique(litter$species)),
    "| species-plot rows:", nrow(litter), "| PCA species:", nrow(traitsonly),
    "\n")

cat("\n== PCA (paper: PC1 48.3%, PC2 38.7%, combined 85.8%) ==\n")
print(round(summary(leaf.PCA)$importance, 3))
print(round(leaf.PCA$rotation, 3))

cat("\n== CTE override (litter$CNratio for CTE, used in models only) ==\n")
print(unique(litter$CNratio[litter$code == "CTE"]))
cat("vs value in PCA:", traitsonly$CNratio[traitsonly$code == "CTE"], "\n")

showCoef <- function(mod, term) {
  ci <- confint(mod, parm = term, method = "Wald")
  cat(sprintf("  %-8s %-8s est = %8.4f  95%% CI = [%.4f, %.4f]\n",
              deparse(substitute(mod)), term, coef(summary(mod))[term, 1],
              ci[1], ci[2]))
}

cat("\n== Trait -> SIV component models (sqrt(iv.wet) ~ trait) ==\n")
cat("Paper: SLA -0.29 [-0.41,-0.17]; PC1 0.07 [0.04,0.11];",
    "LT 1.20 [0.44,1.95]; C:N 0.009 [0.004,0.013]\n")
showCoef(D1c, "SLA")
showCoef(D5c, "PC1")
showCoef(D2c, "LT")
showCoef(D4c, "CNratio")
showCoef(D3c, "LDMC")
showCoef(D6c, "PC2")
cat("R2 (%) [paper: SLA 24.10, PC1 21.82, LT 12.49, CN 18.25, PC2 6.06]:\n")
print(round(100 * sapply(list(SLA = D1c, PC1 = D5c, LT = D2c, CN = D4c,
                              PC2 = D6c),
                         \(m) summary(m)$r.squared), 2))
cat("SIV -> LDMC R2 (%) [paper 3.77]:",
    round(100 * summary(S3c)$r.squared, 2), "\n")

cat("\n== Canopy production model L7: loglitter ~ logba + iv.wet ==\n")
cat("Paper: SIV -1.67 [-4.37, 1.04]; BA 0.93 [0.61, 1.25]\n")
showCoef(L7, "iv.wet")
showCoef(L7, "logba")
print(VarCorr(L7))
