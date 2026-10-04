# Stage 1/2 smoke test: run the historical main script VERBATIM (lines 4-513
# of scripts/historical/, = commit 74befc2) through SEM model selection with
# piecewiseSEM 2.1.2, and compare with Table 1 and Fig. 5 of the paper.
#
# piecewiseSEM 2.1.2 (CRAN, 2020-12-09) is installed in the project-local
# library .Rlib/ (gitignored); see AGENTS.md.
# Run from the repository root:  Rscript audit/05_smoke_test_sem.R

.libPaths(c(normalizePath(".Rlib"), .libPaths()))
stopifnot(packageVersion("piecewiseSEM") == "2.1.2")

script <- readLines(
  "scripts/historical/leaf traits updated with spp by plot and iv wet.R")
run <- script[4:513]
run[91 - 3] <- "# pairs.cor() skipped"  # used before definition; plot only
run <- sub('read.csv("', 'read.csv("data/', run, fixed = TRUE)

pdf(NULL)
failed <- character()
sink(tempfile())  # silence the many summary() prints
for (e in parse(text = run)) {
  res <- try(eval(e, envir = globalenv()), silent = TRUE)
  if (inherits(res, "try-error"))
    failed <- c(failed, substr(deparse(e)[1], 1, 60))
}
sink()
cat("Expressions that failed (expected: SSI.stem models only):\n")
print(failed)

cat("\nModels:", length(mods), "| with AICc:", sum(!is.na(AIC.list)), "\n")

cat("\n== Table 1 reproduction: top 12 ranks ==\n")
top <- final.table[1:12, c("Biol.Hypothesis", "R2", "Fisher.C", "AICc",
                           "dAICc", "weights")]
top$R2 <- gsub("\n", " / ", top$R2)
top$Fisher.C <- round(top$Fisher.C, 3)
top$formula <- sapply(final.table$formula[1:12], function(f)
  gsub("\\s+", " ", paste(f, collapse = " ")))
print(top, right = FALSE)

cat("\nRank of 1st plot-type model (paper: rank 25, dAICc 11.51, w 0.001):\n")
pt <- grep("plotwet", final.table$formula)[1]
print(final.table[pt, c("Biol.Hypothesis", "AICc", "dAICc", "weights")])
cat("  rank", pt, "\n")

cat("\n== Fig. 5 summed weights (paper: vi 61.3, i 19.2, vii 9.3, ii 7.8) ==\n")
print(round(100 * sort(SoW, decreasing = TRUE), 1))

cat("\n== Models whose summary() fails, by rank ==\n")
for (r in seq_along(ranked.index)) {
  i <- ranked.index[r]
  s <- try(summary(mods[[i]], .progressBar = FALSE), silent = TRUE)
  if (inherits(s, "try-error"))
    cat(sprintf("rank %3d (mods[[%3d]], BH %s, dAICc %.2f): %s", r, i, BH[i],
                dAIC[i], conditionMessage(attr(s, "condition"))), "\n")
}

cat("\n== Ranks 11-30 from component formulae (summary()-independent) ==\n")
for (r in 11:30) {
  i <- ranked.index[r]
  f <- sapply(mods[[i]][sapply(mods[[i]], inherits, c("lm", "merMod"))],
              \(m) deparse1(formula(m)))
  f <- gsub(" + (1 | plot) + (logba | species)", "", f, fixed = TRUE)
  cat(sprintf("%3d mods[[%3d]] BH %-4s AICc %.2f d %.2f w %.3f | %s\n", r, i,
              BH[i], AIC.list[i], dAIC[i], unranked.weights[i],
              paste(f, collapse = "; ")))
}
