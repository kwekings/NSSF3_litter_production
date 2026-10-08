# Download the NEA daily weather files that scripts/historical/GAMM only.R
# reads (l. 33-63): S122 Khatib (temperature), S69 Upper Peirce and S40 Mandai
# (rainfall), monthly DAILYDATA_<station>_<YYYYMM>.csv files from
# weather.gov.sg (historical daily records). Same URL pattern as
# TFERP/LTFEM_aux/docs/ElNino/ElNino.Rmd. The files are stored in that sister
# repo, one folder per station (owner's decision, 2026-10-06; committed there
# in b4b9af8). Files already on disk are skipped.
# Run from the repository root:  Rscript audit/07_download_nea.R

dir_nea <- "../../TFERP/LTFEM_aux/data/nea/weather.gov.sg/"
stopifnot(dir.exists(dir_nea))
stations <- c(S122 = "S122_khatib", S69 = "S69_upper_peirce_reservoir",
              S40 = "S40_mandai")

# Which months the historical code reads: months[1:n], where
# months <- c("201806", "201807", unique(tll$date2)) and n = number of unique
# date2 (the vector is offset by two). Download the whole span to be safe.
tll <- read.csv("data/data for GAMM.csv")
months_hist <- c("201806", "201807", gsub("-", "", unique(tll$date2)))
months_hist <- months_hist[seq_along(unique(tll$date2))]
cat("Months read by the historical code:\n")
print(months_hist)

span <- format(seq(as.Date("2018-06-01"), as.Date("2020-08-01"),
                   by = "month"), "%Y%m")
for (stn in names(stations)) {
  dir_stn <- file.path(dir_nea, stations[[stn]])
  dir.create(dir_stn, showWarnings = FALSE)
  for (ym in span) {
    dest <- file.path(dir_stn, sprintf("DAILYDATA_%s_%s.csv", stn, ym))
    if (file.exists(dest)) next
    url <- paste0("https://www.weather.gov.sg/files/dailydata/",
                  basename(dest))
    ok <- !inherits(try(suppressWarnings(
      download.file(url, dest, mode = "wb", quiet = TRUE)),
      silent = TRUE), "try-error")
    if (!ok) unlink(dest)
    if (!ok) message(stn, " ", ym, ": not available")
  }
}
cat("\nFiles on disk for", span[1], "to", tail(span, 1), "by station:\n")
print(sapply(stations, function(s)
  sum(sub("^DAILYDATA_.*_(\\d{6})\\.csv$", "\\1",
          list.files(file.path(dir_nea, s))) %in% span)))
