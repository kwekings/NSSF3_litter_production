# Stage 2 smoke test: run the historical GAMM code VERBATIM (from "START HERE"
# in scripts/historical/GAMM only.R, = commit 74befc2) on the committed
# data/data for GAMM.csv, and compare with the published Fig. 3 results:
#   wet/dry multiplier x2.01 [1.57, 2.57]; temperature 0.40 [0.27, 0.53];
#   adj R2 ~ 75%; best model dAICc 5.47 ahead, weight > 90%.
#
# gamm4 0.2-6 (the version cited in the paper) is in the project-local
# library .Rlib/ (gitignored); see AGENTS.md.
# Run from the repository root:  Rscript audit/06_gamm.R

.libPaths(c(normalizePath(".Rlib"), .libPaths()))
cat("gamm4", format(packageVersion("gamm4")),
    "| MuMIn", format(packageVersion("MuMIn")),
    "| mgcv", format(packageVersion("mgcv")),
    "| lme4", format(packageVersion("lme4")), "\n")

script <- readLines("scripts/historical/GAMM only.R")
start <- grep("### START HERE ###", script, fixed = TRUE)
# Packages (l. 6-8), then the model-selection and best-model code up to
# summary(best$mer); plots are discarded. setwd() is replaced by the data path.
run <- c(script[6:8], script[(start + 1):(start + 44)])
run <- sub("^setwd\\(.*$", "# setwd() skipped", run)
run <- sub('read.csv("data for GAMM.csv"',
           'read.csv("data/data for GAMM.csv"', run, fixed = TRUE)

pdf(NULL)  # the historical code calls plot(); discard
failed <- character()
sink(tempfile())  # silence the summary() prints
for (e in parse(text = run)) {
  res <- try(print(eval(e, envir = globalenv())), silent = TRUE)
  if (inherits(res, "try-error"))
    failed <- c(failed, substr(deparse(e)[1], 1, 60))
}
sink()
cat("Expressions that failed:\n")
print(failed)

cat("\n== Dredge: models within 10 AICc ==\n")
print(dredged[dredged$delta < 10, ])

cat("\n== Best model (paper: Condition + temp) ==\n")
sg <- summary(best$gam)
print(sg$p.table)
cat("adj R2 (gam):", round(sg$r.sq, 4), "\n")

# Effects on the log scale; the paper reports the wet/dry multiplier on the
# response scale and the temperature coefficient as is. Wald 95% CIs.
est <- sg$p.table[, "Estimate"]
se <- sg$p.table[, "Std. Error"]
ci <- cbind(est, lwr = est - 1.96 * se, upr = est + 1.96 * se)
# The reference level is D (dry), so the paper's x2.01 is dry relative to
# wet, i.e. exp(-ConditionW).
cond <- exp(-ci["ConditionW", c("est", "upr", "lwr")])
names(cond) <- c("est", "lwr", "upr")
cat("\nDry vs wet multiplier, exp(-ConditionW), Wald",
    "(paper x2.01 [1.57, 2.57]):\n")
print(round(cond, 3))
cat("temp, Wald (paper 0.40 [0.27, 0.53]):\n")
print(round(ci["temp", ], 3))
cat("Profile CIs on best$mer (log scale):\n")
print(round(confint(best$mer, method = "profile",
                    parm = c("XConditionW", "Xtemp")), 3))

# Historical l. 215: r.squaredLR(best) errors on the gamm4 list under the
# current MuMIn. l. 216 uses the dredge's top model (a uGamm fit), which
# works; its adjusted value is the paper's "adj R2 ~ 75%".
cat("\nr.squaredLR(top dredge model) (paper adj R2 ~ 75%):\n")
print(MuMIn::r.squaredLR(get.models(dredged, subset = 1)[[1]]))
