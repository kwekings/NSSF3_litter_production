# Historical scripts (do not edit)

These are byte-for-byte copies of the scripts as first committed to this
repository in commit `74befc2` (2025-07-31). They are the closest available
record of the code that produced the results in Lam et al. (2022,
*Ecosystems* 25: 1006–1019). They are kept unmodified as the reference for
reproducing the published analysis. The copies in `scripts/` were edited
after 2025-07-31 and no longer represent the published analysis.

The scripts were run interactively and do not run from top to bottom as-is.
They also expect a flat working directory (`setwd()` to a private folder) and
some inputs from outside this repository. Run order:

1. `part I (updated with spp by plot).R` writes
   `litter.production w species-plot CAA Feb21.csv`.
2. `leaf traits updated with spp by plot and iv wet.R` produces the leaf
   trait PCA, the component models, the structural equation models and the
   model selection.
3. `GAMM only.R` produces the phenology GAMM. It uses objects created by
   part I and daily weather files that are not in this repository.

Known issue: the species-level leaf C:N ratio is computed in part I as the
mean over all rows of `CHNS v3.csv`, which mixes fresh and senesced leaf
samples. See the root `README.md`.
