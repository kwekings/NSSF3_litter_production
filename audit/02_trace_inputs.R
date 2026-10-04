# Stage 1 audit: trace the other historical inputs against the intermediate
# files. Read-only with respect to data/; prints checks only.
# Run from the repository root:  Rscript audit/02_trace_inputs.R

library(tidyverse)

feb21 <- read.csv("data/litter.production w species-plot CAA Feb21.csv")
traitSIV <- read.csv("data/leaf traits and SIV.csv")

# 1. Species trait means, rebuilt as in the main script at 74befc2 ----
HJtraits <- read.csv("data/HJ leaf data full.csv") |>
  mutate(across(c(thickness, ldmc, SLA), \(x) suppressWarnings(as.numeric(x))),
         code = toupper(substr(indiv, 1, 3)))
PJtraits <- read.csv("data/PJ leaf dat.csv", check.names = FALSE)
# historical: rowwise mean of the three thickness readings, then species mean
PJtraits$meanlt <- rowMeans(PJtraits[, 6:8], na.rm = TRUE)

traits_HJ <- HJtraits |>
  group_by(code) |>
  summarise(SLA = mean(SLA, na.rm = TRUE),
            LDMC = mean(ldmc, na.rm = TRUE) / 1000,
            LT = mean(thickness, na.rm = TRUE),
            n_leaf = n())
traits_PJ <- PJtraits |>
  group_by(code = SPECIES) |>
  summarise(SLA = mean(`SLA_cm2.g-1`, na.rm = TRUE),
            LDMC = mean(LDMC, na.rm = TRUE),
            LT = mean(meanlt),
            n_leaf = n())
cat("\nCodes in both HJ and PJ (would duplicate rows in rbind/merge):\n")
print(intersect(traits_HJ$code, traits_PJ$code))

traits <- bind_rows(HJ = traits_HJ, PJ = traits_PJ, .id = "source")
cat("\nTrait means vs data/leaf traits and SIV.csv (max abs diff):\n")
traits |>
  inner_join(traitSIV, by = "code", suffix = c("", "_tab")) |>
  summarise(n = n(),
            SLA = max(abs(SLA - SLA_tab)),
            LDMC = max(abs(LDMC - LDMC_tab)),
            LT = max(abs(LT - LT_tab))) |>
  print()
cat("Codes with traits but not in leaf traits and SIV.csv: ",
    setdiff(traits$code, traitSIV$code), "\n")
cat("Trait sources:\n")
print(traits |> filter(code %in% traitSIV$code) |> count(source))

# 2. SIV / SSI from ecophysio_traits/raw_data/SSI Jan21.csv ----
ssiPath <- "../ecophysio_traits/raw_data/SSI Jan21.csv"
if (file.exists(ssiPath)) {
  SSI <- read.csv(ssiPath)
  chk <- feb21 |>
    distinct(species, SSI, iv.wet) |>
    left_join(SSI, by = c("species" = "X"), suffix = c("", "_jan21"))
  cat("\nSpecies in Feb21 not matched in SSI Jan21.csv: ",
      chk$species[is.na(chk$iv.wet_jan21)], "\n")
  cat("Max abs diff SSI vs ssi.ba:",
      max(abs(chk$SSI - chk$ssi.ba), na.rm = TRUE),
      "; iv.wet:", max(abs(chk$iv.wet - chk$iv.wet_jan21), na.rm = TRUE), "\n")
  cat("Madhuca sp. in SSI Jan21: ")
  print(SSI[SSI$X == "Madhuca sp.", ])
  cat("Columns:", names(SSI), "\n")
}
cat("SIV in leaf traits table == iv.wet in Feb21? max abs diff:",
    feb21 |>
      distinct(code, iv.wet) |>
      inner_join(traitSIV, by = "code") |>
      with(max(abs(iv.wet - SIV))), "\n")

# 3. Litter collection, start and duration from the lamina data ----
lamina <- read.csv("data/Cleaned lamina+PJ.csv") |>
  mutate(date1 = as.POSIXct(General.date, format = "%d/%m/%Y", tz = "GMT"))
coll <- lamina |>
  group_by(plot = Plot, species) |>
  summarise(litter.collection_re = sum(Dry.Mass), .groups = "drop")
dur <- lamina |>
  group_by(species) |>
  summarise(start_re = as.character(as.Date(min(date1, na.rm = TRUE))),
            dmax = max(date1, na.rm = TRUE), dmin = min(date1, na.rm = TRUE)) |>
  mutate(duration_re = as.numeric(
    ifelse(dmax == dmin,
           difftime(max(lamina$date1, na.rm = TRUE), dmin, units = "days"),
           difftime(dmax, dmin, units = "days"))))
chk <- feb21 |>
  left_join(coll, by = c("plot", "species")) |>
  left_join(dur, by = "species")
cat("\nLitter collection max abs diff:",
    max(abs(chk$litter.collection - chk$litter.collection_re)),
    "; duration max abs diff:", max(abs(chk$duration - chk$duration_re)),
    "; start mismatches:", sum(chk$start != chk$start_re), "\n")

# 4. Basal area: is ten.ba recoverable from NSSF2_main tree data? ----
treePath <- "../NSSF2_main/Data/NSSF2trees_160324.csv"
if (file.exists(treePath)) {
  tree <- read.csv(treePath) |>
    mutate(plot = paste0("Q", plot),
           ba = pi * (suppressWarnings(as.numeric(dbh)) / 2)^2) |>
    group_by(plot, species) |>
    summarise(ba_nssf2 = sum(ba, na.rm = TRUE), .groups = "drop")
  chk <- feb21 |> left_join(tree, by = c("plot", "species"))
  cat("\nten.ba vs NSSF2trees_160324 dbh: exact matches",
      sum(abs(chk$ten.ba - chk$ba_nssf2) < 1e-6, na.rm = TRUE),
      "of", nrow(chk), "; unmatched", sum(is.na(chk$ba_nssf2)), "\n")
  print(head(select(chk, plot, species, ten.ba, ba_nssf2), 8))
}

# 5. HJ leaf data full.csv vs ecophysio_traits 'leaf soft traits.csv' ----
softPath <- "../ecophysio_traits/raw_data/leaf soft traits.csv"
if (file.exists(softPath)) {
  soft <- read.csv(softPath)
  cat("\nRows: HJ full", nrow(HJtraits), "; leaf soft traits", nrow(soft),
      "\n")
  cmp <- HJtraits |>
    select(indiv, area, fresh_mass, dry_mass, thickness, ldmc, SLA) |>
    left_join(soft, by = c("indiv" = "twig"), suffix = c("", "_soft"))
  for (v in c("area", "fresh_mass", "dry_mass", "thickness", "ldmc", "SLA")) {
    d <- abs(as.numeric(cmp[[v]]) - as.numeric(cmp[[paste0(v, "_soft")]]))
    cat(sprintf("  %-10s rows differing (>1e-6): %d\n", v,
                sum(d > 1e-6, na.rm = TRUE)))
  }
}
