# Stage 1 audit: can the historical species-level C:N values be reconstructed
# from the fresh and senesced CHNS rows?
#
# Read-only with respect to data/; writes only to audit/.
# Run from the repository root:  Rscript audit/01_cn_diagnostic.R

library(tidyverse)

# 1. Check the 2025 split against the original CHNS v3.csv ----
# The original file is no longer in the working tree; read it from the commit
# that added it (71eb215).

chns_orig <- system('git show "71eb215:data/CHNS v3.csv"', intern = TRUE) |>
  paste(collapse = "\n") |>
  I() |>
  read.csv(text = _)

chns_fresh <- read.csv("data/CHNS v3_fresh.csv")
chns_senesced <- read.csv("data/CHNS v3_senesced.csv")

# Spreadsheet rows 2-74 = data rows 1-73; rows 75-126 = data rows 74-125
stopifnot(
  nrow(chns_orig) == 125,
  isTRUE(all.equal(chns_orig[1:73, ], chns_senesced, check.attributes = FALSE)),
  isTRUE(all.equal(chns_orig[74:125, ], chns_fresh,
                   check.attributes = FALSE))
)
message("Split check passed: senesced = rows 2-74, fresh = rows 75-126.")

# 2. Species-level means under each derivation ----
# The historical code (part I, line 226 at commit 74befc2) was
#   cn <- with(chns, tapply(Ratio, Species, mean))
# i.e. the mean of per-sample ratios over ALL rows, fresh and senesced alike.

cn <- bind_rows(fresh = chns_fresh, senesced = chns_senesced,
                .id = "type_leaf") |>
  group_by(Species) |>
  summarise(n_fresh = sum(type_leaf == "fresh"),
            n_senesced = sum(type_leaf == "senesced"),
            cn_fresh = mean(Ratio[type_leaf == "fresh"]),
            cn_senesced = mean(Ratio[type_leaf == "senesced"]),
            cn_combined = mean(Ratio)) |>
  mutate(across(starts_with("cn_"), \(x) ifelse(is.nan(x), NA, x)))

# 3. Historical values ----
# (a) litter.consol output of part I, carrying CNratio per species
feb21 <- read.csv("data/litter.production w species-plot CAA Feb21.csv") |>
  distinct(species, code, CNratio_feb21 = CNratio)

# (b) species-level table believed to match the published analysis
# (Madhuca tomentosa relabelled Madhuca sp. in this file)
traitSIV <- read.csv("data/leaf traits and SIV.csv") |>
  mutate(species = recode(species, "Madhuca sp." = "Madhuca tomentosa")) |>
  select(species, code_traitSIV = code, CNratio_traitSIV = CNratio)

# 4. Species entering the trait analysis ----
# Main script (74befc2): traits are keyed by 3-letter code from HJ (first three
# letters of `indiv`) and PJ (`SPECIES`); merge() with the C:N codes keeps the
# intersection.
codes_HJ <- read.csv("data/HJ leaf data full.csv")$indiv |>
  substr(1, 3) |>
  toupper() |>
  unique()
codes_PJ <- unique(read.csv("data/PJ leaf dat.csv")$SPECIES)

diag <- cn |>
  full_join(feb21, by = c("Species" = "species")) |>
  full_join(traitSIV, by = c("Species" = "species")) |>
  mutate(in_analysis = code %in% c(codes_HJ, codes_PJ),
         diff_combined_hist = cn_combined - CNratio_feb21,
         diff_fresh_hist = cn_fresh - CNratio_feb21,
         status = case_when(
           is.na(CNratio_feb21) ~ "not in litter data / no C:N",
           n_fresh > 0 & n_senesced > 0 ~ "MIXED fresh+senesced",
           n_fresh > 0 ~ "fresh only",
           n_senesced > 0 ~ "SENESCED only"
         )) |>
  arrange(desc(in_analysis), status, Species)

# The two historical sources must agree where both exist
stopifnot(with(diag, all(abs(CNratio_feb21 - CNratio_traitSIV) < 1e-6,
                         na.rm = TRUE)))
# The historical values must equal the all-rows mean
stopifnot(with(diag, all(abs(diff_combined_hist) < 1e-6, na.rm = TRUE)))
message("Historical CNratio == mean of ALL CHNS rows, for every species.")

write.csv(diag, "audit/cn_diagnostic.csv", row.names = FALSE)

diag |>
  select(Species, code, in_analysis, status, n_fresh, n_senesced,
         cn_fresh, cn_senesced, cn_combined, CNratio_feb21,
         diff_fresh_hist) |>
  mutate(across(where(is.double), \(x) round(x, 2))) |>
  as.data.frame() |>
  print()

# 5. Code collisions in the 3-letter keys used for merging ----
read.csv("data/Cleaned lamina+PJ.csv") |>
  distinct(species) |>
  filter(str_detect(species, "^[A-Z][a-z]+ [a-z-]+")) |>
  mutate(code = toupper(paste0(substr(species, 1, 1),
                               substr(word(species, 2), 1, 2)))) |>
  group_by(code) |>
  filter(n() > 1) |>
  arrange(code) |>
  print()
