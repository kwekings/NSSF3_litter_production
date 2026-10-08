# Rebuild the weather columns of data/data for GAMM.csv (rf, rf_lag, temp,
# temp_lag, temp_max, temp_min) by running the historical weather code
# VERBATIM (scripts/historical/GAMM only.R l. 32-121) on NEA files freshly
# downloaded by audit/07_download_nea.R, and compare with the committed file.
#
# Patches (paths and encoding only):
# - the two D:\ folders become the station folders in the sister repo
#   TFERP/LTFEM_aux/data/nea/weather.gov.sg/;
# - NEA files are read as CP1252, where the missing-value byte 0x97 is an
#   em dash. The historical script holds the same raw 0x97 byte in its
#   missing-value literal, so it is read as CP1252 too, and the historical
#   missing-value handling applies unchanged (missing rain -> 0; missing
#   mean temp -> monthly mean).
# The litter part of tll (Dry.Mass, ba, prod) is taken from the committed
# file, already verified against lamina + trees (AGENTS.md, Sec. 8).
# Run from the repository root:  Rscript audit/08_rebuild_gamm_weather.R

ref <- read.csv("data/data for GAMM.csv")
ref$X <- NULL
tll <- ref[c("Plot", "date1", "date2", "Dry.Mass", "Condition", "ba",
             "prod")]

script <- iconv(readLines("scripts/historical/GAMM only.R"),
                "CP1252", "UTF-8")
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
eval(parse(text = run, encoding = "UTF-8"), envir = globalenv())

cat("Weather days read:", nrow(weather), "|",
    format(min(weather$date1)), "to", format(max(weather$date1)), "\n")
cat("Missing rain set to 0:",
    sum(rainfall$total_rf == "\u2014", na.rm = TRUE), "days\n")

vars <- c("rf", "rf_lag", "temp", "temp_lag", "temp_max", "temp_min")
cmp <- data.frame(var = vars, t(sapply(vars, function(v) {
  d <- tll[[v]] - ref[[v]]
  c(exact = sum(abs(d) < 1e-6, na.rm = TRUE),
    na_rebuilt = sum(is.na(tll[[v]])),
    max_abs_diff = round(max(abs(d), na.rm = TRUE), 3))
})))
cat("\nRebuilt vs committed (", nrow(tll), "rows):\n")
print(cmp, row.names = FALSE)

# Which collection months disagree?
bad <- !Reduce(`&`, lapply(vars, function(v)
  abs(tll[[v]] - ref[[v]]) < 1e-6 & !is.na(tll[[v]])))
cat("\nRows with any mismatch, by collection month:\n")
print(table(substr(tll$date1, 1, 7), bad))
