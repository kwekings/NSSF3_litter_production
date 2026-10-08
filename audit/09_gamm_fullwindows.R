# Sensitivity of the GAMM (Fig. 3) to the truncated weather windows.
#
# The historical weather code (scripts/historical/GAMM only.R l. 32-121) reads
# monthly files via an offset `months` vector that stops at 201906, so the 24
# collections in July 2019 have truncated 30-day windows (down to 1 day).
# Here the same code is run with ONE change: `months` covers 201806-201907,
# so every window is complete. Everything else is historical: Mandai (S40)
# rain for 201903-201906, missing rain -> 0, missing mean temp -> monthly
# mean. Then the historical GAMM code (from "START HERE") is refitted on the
# rebuilt data and compared with the historical fit (audit/06_gamm.R).
#
# Needs the NEA files in TFERP/LTFEM_aux (audit/07_download_nea.R) and
# gamm4 0.2-6 in .Rlib/.
# Run from the repository root:  Rscript audit/09_gamm_fullwindows.R

.libPaths(c(normalizePath(".Rlib"), .libPaths()))

ref <- read.csv("data/data for GAMM.csv")
ref$X <- NULL
tll <- ref[c("Plot", "date1", "date2", "Dry.Mass", "Condition", "ba",
             "prod")]

# The historical script holds the raw CP1252 missing-value byte 0x97, like
# the NEA files; see audit/08_rebuild_gamm_weather.R.
script <- iconv(readLines("scripts/historical/GAMM only.R"),
                "CP1252", "UTF-8")

## Rebuild the weather columns with full windows ----
run <- script[32:121]
# NEA files live in the sister repo TFERP/LTFEM_aux, one folder per station
# (audit/07_download_nea.R). `location` is blanked and the folder goes
# into each file prefix, because the rain loop reads two stations.
dir_nea <- "../../TFERP/LTFEM_aux/data/nea/weather.gov.sg/"
run <- sub('^location <- "D:.*$', 'location <- ""', run)
run <- sub('"DAILYDATA_S122_"',
           'paste0(dir_nea, "S122_khatib/DAILYDATA_S122_")', run,
           fixed = TRUE)
run <- sub('"DAILYDATA_S69_"',
           'paste0(dir_nea, "S69_upper_peirce_reservoir/DAILYDATA_S69_")',
           run, fixed = TRUE)
run <- sub('"DAILYDATA_S40_"',
           'paste0(dir_nea, "S40_mandai/DAILYDATA_S40_")', run,
           fixed = TRUE)
run <- gsub('.csv"))', '.csv"), fileEncoding = "CP1252")', run,
            fixed = TRUE)
# The one substantive change: read every month from 201806 to 201907
i_months <- grep("^months <- ", run)
stopifnot(length(i_months) == 1)
run[i_months] <- paste0('months <- format(seq(as.Date("2018-06-01"), ',
                        'as.Date("2019-07-01"), by = "month"), "%Y%m")')
n_loops <- sum(grepl("for(i in 1:length(unique(tll$date2)))", run,
                     fixed = TRUE))
stopifnot(n_loops == 2)
run <- gsub("for(i in 1:length(unique(tll$date2)))",
            "for(i in seq_along(months))", run, fixed = TRUE)
eval(parse(text = run, encoding = "UTF-8"), envir = globalenv())

cat("Weather days read:", nrow(weather), "|",
    format(min(weather$date1)), "to", format(max(weather$date1)), "\n")

vars <- c("rf", "rf_lag", "temp", "temp_lag", "temp_max", "temp_min")
changed <- Reduce(`|`, lapply(vars, function(v)
  abs(tll[[v]] - ref[[v]]) > 1e-6))
cat("\nRows whose weather values changed:", sum(changed), "\n")
cat("Collection dates affected:\n")
print(table(format(tll$date1[changed])))
cat("\nMean change in the affected rows (rebuilt - historical):\n")
print(round(sapply(vars, function(v)
  mean(tll[[v]][changed] - ref[[v]][changed])), 3))

## Refit the historical GAMM code on the rebuilt data ----
library(mgcv)
library(MuMIn)
library(gamm4)
start <- grep("### START HERE ###", script, fixed = TRUE)
run <- script[(start + 1):(start + 44)]
run <- sub("^setwd\\(.*$", "# setwd() skipped", run)
# Use the rebuilt tll instead of reading the committed file
run <- sub('^tll <- read.csv\\("data for GAMM.csv".*$',
           "# tll rebuilt above", run)
pdf(NULL)  # the historical code calls plot(); discard
failed <- character()
sink(tempfile())  # silence the summary() prints
for (e in parse(text = run)) {
  res <- try(print(eval(e, envir = globalenv())), silent = TRUE)
  if (inherits(res, "try-error"))
    failed <- c(failed, substr(deparse(e)[1], 1, 60))
}
sink()
cat("\nExpressions that failed:\n")
print(failed)

cat("\n== Dredge: models within 10 AICc (full windows) ==\n")
print(dredged[dredged$delta < 10, ])

# The historical code hard-codes `best` as Condition + temp; report it
# whether or not it is still the top-ranked model.
sg <- summary(best$gam)
est <- sg$p.table[, "Estimate"]
se <- sg$p.table[, "Std. Error"]
ci <- cbind(est, lwr = est - 1.96 * se, upr = est + 1.96 * se)
cond <- exp(-ci["ConditionW", c("est", "upr", "lwr")])
r2_top <- MuMIn::r.squaredLR(get.models(dredged, subset = 1)[[1]])

# Historical values from audit/06_gamm.R (same software versions)
comparison <- data.frame(
  quantity = c("top-model weight", "dAICc to 2nd model",
               "dry/wet multiplier", "  lower 95%", "  upper 95%",
               "temp coefficient", "  lower 95%", "  upper 95%",
               "adj R2 (LR, top model)"),
  historical = c(0.908, 5.47, 2.020, 1.575, 2.590, 0.396, 0.264, 0.529,
                 0.752),
  fullWindows = round(c(dredged$weight[1], dredged$delta[2],
                        cond, ci["temp", ],
                        attr(r2_top, "adj.r.squared")), 3),
  paper = c(">0.90", "5.47", "2.01", "1.57", "2.57", "0.40", "0.27",
            "0.53", "~0.75"))
cat("\n== Historical vs full-window fit",
    "(Wald CIs; `best` = Condition + temp) ==\n")
print(comparison, row.names = FALSE)
