# Chain Elongation Validation Framework

Code accompanying:

> **Quantitative Fragility and Ecological Robustness in Machine-Learning-Predicted Chain Elongation: A Multi-Study, Multi-Platform Validation Framework**
> Preetam Banerjee. Water Supply and Bioeconomy Division, Poznań University of Technology. MSCA-LeAD, BIOTWIN/BioRef, WP6 Task 6.5.

## What this is

Microbiome-to-function machine learning models are typically validated within a single experimental system. This project asks what happens when a caproate-prediction model built on six previously published short-read (Illumina) chain-elongation studies plus one newly contributed long-read (Nanopore) dataset is validated properly instead: by study, by individual physical reactor, and across sequencing platforms.

The short version: naive validation reaches R² as high as 0.985; every properly grouped scheme collapses that to at or below zero. The taxa the model identifies as core drivers, however, stay stable across every validation architecture, every algorithm tested, and even across sequencing platforms. That contrast — quantitative fragility alongside ecological robustness — is the paper's central finding.

This repository contains the complete pipeline: raw sequence processing through final figures, with every reported number traceable to a script that produces it.

## Repository structure

```
R/                     Sequence processing, in the order it was actually run
  01a-01f_...           DADA2 processing for each of the six short-read studies
  02_taxonomic_assignment.R
  03_grand_merge.R       Merges the six studies' ASV tables into one matrix
  04_metadata_merge.R    Joins operational/metabolite metadata (Day-0 inoculum
                          duplication logic, NA-vs-zero handling)

python/pipeline/       The validation architectures, in Methods 2.6's order
  01_genus_aggregation.py            ASV -> genus collapse (Methods 2.4)
  02_loso_by_study.py                Leave-one-study-out, all 4 feature tiers
                                      (Table 3b's LOSO row + all of Table 4)
  03_permutation_test_999.py         999-iteration label-permutation control
  04_vessel_grouped_deployability.py GroupKFold by physical reactor (Table 3b)
  05_intrastudy_transferability_table5.py   Final leakage-controlled reactor-pair
                                             transfer (Table 5; Figure 3, stage 3)
  06_cross_platform_transfer.py      Illumina -> Nanopore transfer (Table 3b)
  07a_algorithm_comparison_lightgbm_rf.py        Algorithm selection (Methods 2.7);
  07b_algorithm_comparison_lightgbm_standard.py  the target-leakage demonstration
                                                  behind excluding downstream co-products
  07c_figure3_intermediate_unscaled_xgboost.py   Figure 3, stage 2 (R2=-0.228)
  07d_figure3_stage1_random_forest.py            Figure 3, stage 1 (R2=0.217)

python/figures/        One script per figure, each independently runnable
  figure1_r2_comparison.py           Master validation-collapse bar chart
  figure2_permutation_histogram.py   Loads results/permutation_null_distribution_999.csv
  figure3_duber_progression.py       The three-stage leakage-tightening figure
  figure4_cohort_comparison.py       Illumina vs. Nanopore, real per-sample data
  figure5_shap_genus.py              Real shap.TreeExplainer output, genus-only
  figure_s1_substrate_ladder.py
  figure_s3_intrastudy_transfer.py
  figure_s4_duber_scatter.py          Real recomputation of Table 5's worst fold
  figure_s5_combined_shap.py          Chemistry + genus combined SHAP

results/                Saved output from long-running scripts (see note below)
figures_output/         Where the figure scripts write PNG/PDF output
```

Every figure in the manuscript, main and supplementary, has a corresponding script in `python/figures/`.

## Reproducing this

**Shortcut**: `data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv` is already the final matrix -- if you just want to regenerate figures or rerun validation architectures, skip to step 4.

**Full pipeline from raw sequences:**
1. Install R packages (see `R_packages.txt`) and Python packages (`pip install -r requirements.txt`), using the **exact pinned versions**, not just compatible ones — see the reproducibility note below.
2. Run `R/` in numeric order to go from raw SRA accessions (Table 1) to the merged, pre-genus-collapse matrix. Each script has an "EDIT THIS" marker where you need to set your own local working directory -- these manage multi-stage SRA download and DADA2 output and are not one-command reruns.
3. Run `python/pipeline/01_genus_aggregation.py` (using `data/historical/BIOTWIN_FINAL_ML_MATRIX.csv` and `ASV_Taxonomy_Master.csv`) to reproduce `data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv`, the 184-sample x 170-genus matrix every other script reads.
4. Run the remaining `python/pipeline/` scripts in numeric order to reproduce Tables 3a, 3b, 4, and 5. Scripts 07a-07d use the historical snapshots in `data/historical/` instead (see each script's header).
5. Run any `python/figures/` script independently to regenerate that figure -- run from inside `python/figures/` or `python/pipeline/`, since paths are relative to each script's own location.

`python/pipeline/03_permutation_test_999.py` takes roughly 35 minutes (999 iterations x 7-fold LOSO x XGBoost fit) and saves its full result array to `results/permutation_null_distribution_999.csv`. `figure2_permutation_histogram.py` reads from that file rather than recomputing — run the permutation script first.

## A note on exact reproducibility

XGBoost does not guarantee bit-identical output across library versions, even with a fixed `random_state`; internal tie-breaking and numerics can differ between releases. Every number in this repository was generated with the exact package versions pinned in `requirements.txt` and `R_packages.txt`. Re-running with different versions will very likely reproduce the same qualitative pattern (naive validation inflated, grouped validation collapsed, core taxa stable) but may not reproduce every reported decimal exactly — this was confirmed directly during development (see commit history / project log) and is a known characteristic of multi-threaded gradient boosting, not a bug in this code.

## Data availability

Raw sequencing data for the six short-read studies are available under their original NCBI/ENA accessions (Table 1 of the manuscript). The Nanopore dataset (Prusak et al., 2026) is available from the corresponding collaborator upon reasonable request pending its own independent publication.

## Data

```
data/
  BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv   The canonical 184-sample matrix.
                                              Used by every script except
                                              07a/07b/07c/07d below.
  historical/
    BIOTWIN_GENUS_ML_MATRIX.csv              Used by 07c and 07d -- an
                                              intermediate snapshot predating
                                              corrections folded into the
                                              final matrix. Included for
                                              provenance and to keep Figure 3's
                                              full three-stage progression
                                              reproducible, not as a canonical
                                              dataset for any main result.
    BIOTWIN_PRUNED_ML_MATRIX.csv             Used by 07a/07b (algorithm
                                              selection, ASV-level, predates
                                              genus aggregation).
    BIOTWIN_FINAL_ML_MATRIX.csv,
    ASV_Taxonomy_Master.csv                  Inputs to 01_genus_aggregation.py.
```

All four dataset snapshots are provided so the complete history behind Figure 3 -- not just its final, headline number -- stays runnable, in keeping with this repository's whole premise: every reported number traceable to a script that produces it.

## Resolved during development (kept here for the record)

- Figure 3's full three-stage progression (R²=0.217 → −0.228 → −5.869) is now backed end to end: `07d_figure3_stage1_random_forest.py`, `07c_figure3_intermediate_unscaled_xgboost.py`, and `05_intrastudy_transferability_table5.py` respectively.
- The manuscript's description of the middle stage was corrected from "XGBoost without downstream co-products" to "an earlier XGBoost configuration lacking both explicit downstream-co-product exclusion and confirmed training-only scaling" -- the original code retains Valerate and Heptanoate as usable features, visible directly in that script's own SHAP output.

## License

MIT (see `LICENSE`). Cite the manuscript above if you use this code or reuse this validation framework.
