# Download the NEA daily weather files that scripts/historical/GAMM only.R
# reads (l. 33-63): S122 Khatib (temperature), S69 Upper Peirce and S40 Mandai
# (rainfall), monthly DAILYDATA_<station>_<YYYYMM>.csv files from
# weather.gov.sg (historical daily records). Same URL pattern as
# TFERP/LTFEM_aux/docs/ElNino/ElNino.Rmd. Files already on disk are skipped.
# Run from the repository root:  Rscript audit/07_download_nea.R

dir_nea <- "data/nea"
dir.create(dir_nea, showWarnings = FALSE)

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
for (stn in c("S122", "S69", "S40")) {
  for (ym in span) {
    dest <- file.path(dir_nea, sprintf("DAILYDATA_%s_%s.csv", stn, ym))
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
cat("\nFiles on disk:", length(list.files(dir_nea)), "\n")
