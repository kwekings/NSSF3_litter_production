This repository contains udpated data and scripts from the publication:

Lam, W.N., P.J. Chan, Y.Y. Ting, H.J. Sim, J.J. Lian, R. Chong, Nur Estya Rahman, L.W.A. Tan, Q.Y. Ho, Z. Chiam, Srishti Arora, Sorain J. Ramchunder, K.S.-H. Peh, Y. Cai, K.Y. Chong. (2022) Habitat specialization mediates the influence of leaf traits on canopy productivity. Ecosystems 25: 1006–1019. [[website](https://onlinelibrary.wiley.com/share/author/9JQP3IWSSCNUYDRICHSQ?target=10.1111/btp.12913)] [[pdf](https://www.dropbox.com/scl/fi/pfzt1lozsbufmu32nw3q6/lam-et-al-2022-habitat-adaptation-mediates-the-influence-of-leaf-traits-on-canopy-productivity-from-a-tropical-freshwater-swamp-forest.pdf?rlkey=7bfvysjelq4i5mbd41xntf5dm&dl=0)]

The data was deposited on [Figshare](https://figshare.com/articles/dataset/Habitat_adaptation_mediates_the_influence_of_leaf_traits_on_canopy_productivity_evidence_from_a_tropical_freshwater_swamp_forest/15163782?file=29128737) in accordance with the journal's data archival policy. Here, we also provide the scripts that have been further edited since the publication of the paper.

Changes to the data from the Figshare version:

- _Madhuca tomentosa_ has been corrected to _Madhuca_ sp.

Leaf lamina C and N content are divided into:
* Fresh leaf C and N, of which the data from 22 species are shared with the same species in Lam et al. New Phytologist (see [GitHub repository](https://github.com/wengngai/ecophysio_traits) “[CN ratio.csv](https://github.com/wengngai/ecophysio_traits/blob/main/raw_data/CN%20ratio.csv)”) but
  * 16 species are in the Lam et al. New Phytologist paper’s data but not the Ecosystems paper’s data.
  * 4 species in the Ecosystems paper’s data not in the Lam et al. New Phytologist paper’s data.
* Senesced leaf C and N. There are at two rows for most species because we sent in two samples for analyses, except:
  * Adenanthera malayana where a third was sent in because one of the two samples differed more than expected.
  * There were multiple _Madhuca_ sp. because some were originally thought to be _Gluta wallichii_ but the identification was corrected later.

## Historical scripts

`scripts/historical/` holds byte-for-byte copies of the scripts as first committed to this repository (commit `74befc2`, 31 July 2025). They are the closest available record of the code that produced the published results, and are kept unmodified as the reference for reproducing the published analysis. Do not edit them. The copies directly under `scripts/` were edited after that date and no longer represent the published analysis.

The historical scripts were run interactively and do not run from top to bottom as-is; they also expect a flat working directory (`setwd()` to a private folder). Run order:

1. `part I (updated with spp by plot).R` writes `litter.production w species-plot CAA Feb21.csv`.
2. `leaf traits updated with spp by plot and iv wet.R` produces the leaf trait PCA, the component models, the structural equation models and the model selection.
3. `GAMM only.R` produces the phenology GAMM, reading `data for GAMM.csv`.

Known issue: in `part I`, the species-level leaf C:N ratio is computed as the mean over all rows of `CHNS v3.csv`, which mixes fresh and senesced leaf samples (now split into `data/CHNS v3_fresh.csv` and `data/CHNS v3_senesced.csv`).

(To be continued...)
