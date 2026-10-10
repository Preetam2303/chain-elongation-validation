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
  02_grand_merge.R       Merges the six studies' ASV tables into one matrix and
                          removes chimeras
  03_taxonomic_assignment.R  Assigns SILVA v138.1 taxonomy to the merged ASVs
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

results_reference.json  Every reported value, the script that produces it, and
                        how exactly it should reproduce. Start here.
```

Every figure in the manuscript, main and supplementary, has a corresponding script in `python/figures/`.

## Reproducing this

**Shortcut**: `data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv` is the final matrix -- if you just want to regenerate figures or rerun validation architectures, skip to step 4.

**Full pipeline from raw sequences:**
1. Install R packages (see `R_packages.txt`) and Python packages (`pip install -r requirements.txt`), using the **exact pinned versions**, not just compatible ones — see the reproducibility note below.
2. Run `R/` in numeric order to go from raw SRA accessions (Table 1) to the merged, pre-genus-collapse matrix. Each script has an "EDIT THIS" marker where you need to set your own local working directory -- these manage multi-stage SRA download and DADA2 output and are not one-command reruns.
   For the primer-trimmed rerun, `python/pipeline/00_prevalence_filter_genus.py` then applies the 5% prevalence filter (7 of 130 samples) and sums the kept ASVs by genus, writing `BIOTWIN_PRUNED_ML_MATRIX.csv` and `BIOTWIN_GENUS_ML_MATRIX.csv` into the dated `primer_trimmed_*` run folder.
3. Run `python/pipeline/01_genus_aggregation.py` (using `data/historical/BIOTWIN_FINAL_ML_MATRIX.csv` and `ASV_Taxonomy_Master.csv`) to reproduce `data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv`, the 184-sample x 170-genus matrix every other script reads.
4. Run the remaining `python/pipeline/` scripts in numeric order to reproduce Tables 3a, 3b, 4, and 5. Scripts 07a-07d use the historical snapshots in `data/historical/` instead (see each script's header).
5. Run any `python/figures/` script independently to regenerate that figure -- run from inside `python/figures/` or `python/pipeline/`, since paths are relative to each script's own location.
6. Check what you get against `results_reference.json`, which lists every reported value alongside the tolerance you should expect it to reproduce within.

`python/pipeline/03_permutation_test_999.py` takes roughly 35 minutes (999 iterations x 7-fold LOSO x XGBoost fit) and saves its full result array to `results/permutation_null_distribution_999.csv`. `figure2_permutation_histogram.py` reads from that file rather than recomputing — run the permutation script first.

## A note on exact reproducibility

**Read this before reporting a mismatch.** Which numbers reproduce exactly is not uniform across this repository, and the split is predictable.

Every model in this project uses `subsample=0.8, colsample_bytree=0.8`. That stochastic row/column subsampling draws from a pseudo-random stream whose *consumption pattern* is not guaranteed stable across XGBoost builds or platforms. Fixing `random_state=42` therefore gives bit-identical results **within** a fixed environment, but not necessarily **across** environments. This is not thread nondeterminism: results were verified identical across `n_jobs` values of 1, 2, 4 and -1, and identical across repeated runs.

Independent re-execution on a different OS and Python build (Linux, Python 3.12) with the pinned package versions gave:

| Reproduces exactly (no subsampling in the fitted model) | Environment-sensitive (subsampling active) |
|---|---|
| Data construction end to end — 18,870 → 4,027 ASVs, 65 genera, 124 Nanopore genera, 184 × 170 final matrix, all metabolite conversions | Table 3a naive baselines |
| LOSO ablation ladder, all four tiers (Table 4) | GroupKFold by vessel, both modes |
| Both cross-platform transfers (Table 3b) | LOSO SHAP top-25 variant |
| Label-permutation observed value (R² = −0.095) | Reactor-pair transfers (Table 5, Figure 3) |
| Table 1 sample counts, Table 2 prevalence, replicate structure | |
| Figure 3 stage 1 (sklearn RandomForest) | |

The permutation test's *p*-value is a Monte-Carlo estimate over 999 shuffles with a standard error of roughly 0.008; it is seeded per iteration (`random_state=iteration`) and reproduces exactly under that seeding, but is not a fixed quantity under a different shuffle sequence.

Subsampling was retained rather than removed because it is standard regularization and materially improves the fitted models. The cost is the environment sensitivity documented above.

`results_reference.json` in the repository root records every reported value together with its tolerance class, so a mismatch can be checked against expectation rather than guessed at. Values marked `exact` should match to three decimals anywhere; values marked `environment_sensitive` should match in sign and conclusion only.

**Environment used for every reported number:**

```
Python 3.13.9
Windows
package versions exactly as pinned in requirements.txt and R_packages.txt
```

## Data availability

Raw sequencing data for the six short-read studies are available under their original NCBI/ENA accessions (Table 1 of the manuscript). The Nanopore dataset (Prusak et al., 2026) is available from the corresponding collaborator upon reasonable request pending its own independent publication.

## Data

```
data/
  BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv   The canonical 184-sample matrix.
                                              Used by every script except
                                              07a/07b/07c/07d below.
  study_provenance.csv                       One row per short-read study: Paper_ID,
                                              R script, paper, sequence accession used
                                              and the one the paper states, sample counts.
  corrections/
    excluded_samples.csv                     Samples removed from the modelling data,
                                              each with its reason and source.
    metadata_corrections.csv                 Cell-level fixes (old value, new value,
                                              source, reason), applied by
                                              R/04_metadata_merge.R (Illumina) and at
                                              the Nanopore merge (Hanna_2025 rows).
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

## Primer-trimmed rerun: chimera removal

After primer trimming, chimera removal (`removeBimeraDenovo`, consensus method, DADA2's default `minFoldParentOverAbundance = 1.5`) still removes far more reads in three studies than in the other three. The setting was kept at the default. The checks behind that choice are in `R/diagnostics/`, and their outputs are in `results/chimera_check_2026-10-08/`.

- `chimera_check_per_study.R`: running consensus once per study removes the same share of reads as running it once on the merged table (for example, 42.2% vs 42.0% for Brodowski 2022), so merging the studies is not the cause.
- `chimera_fold_sensitivity.R`: percent of reads removed per study when a chimera's parents must be at least 1.5, 4 or 8 times as abundant as it.

| Study | 1.5× (used) | 4× | 8× |
|---|---|---|---|
| Brodowski_PRJNA715197 (external acetate) | 39.6 | 32.4 | 26.4 |
| Brodowski_2022 | 42.2 | 37.7 | 32.2 |
| Duber_2025 | 37.0 | 34.8 | 28.5 |
| Brodowski_2025 | 16.1 | 15.9 | 15.7 |
| Duber_2022 | 12.0 | 10.7 | 9.9 |
| Duber_2020 | 16.4 | 15.8 | 14.2 |

Even at 8× the first three studies lose 26-32% of reads, so most of their loss is not an artefact of the parent-abundance threshold. In one Brodowski 2022 sample (SRR18962128), the most abundant flagged sequences were 1-2 base variants of abundant Enterobacteriaceae sequences. BLAST could not separate them from their parents: flagged sequences and parents all matched type strains at 99.3-99.8% identity over their full length.

`chimera_fold_genus_check.R` then compared the genus picture of those three studies at 1.5× and 8× (`top_genera_fold_1.5_vs_8.csv`, mean per-sample relative abundance). The genus picture does not change. In each study the top two genera keep their order, no named top-10 genus moves more than two places, and none changes by more than 4 percentage points. Caproiciproducens stays first or second in all three studies: 29.9% vs 31.7% (external acetate), 28.8% vs 28.9% (Brodowski 2022) and 28.2% vs 29.4% (Duber 2025). The largest shifts are Acinetobacter in the external-acetate study (36.4% to 32.7%) and the unclassified bins, which gain 0.9-1.8 points at 8×.

## Resolved during development (kept here for the record)

- Figure 3's full three-stage progression (R²=0.217 → −0.228 → −5.869) is backed end to end: `07d_figure3_stage1_random_forest.py`, `07c_figure3_intermediate_unscaled_xgboost.py`, and `05_intrastudy_transferability_table5.py` respectively.
- The manuscript's description of the middle stage has been corrected. It previously read "XGBoost without downstream co-products", which the code contradicts: `07c` drops only sparse/junk columns and leaves Valerate, Isobutyrate and Heptanoate available as predictors — all three appear in that script's own SHAP output, alongside DAY as the single dominant feature. Section 3.5 now describes the stage accurately.
- Table 3b previously reported the headline LOSO result as R² = −0.095 with no standard deviation, while Table 4 reported the same tier-4 value as −0.095 ± 0.387. Table 3b, Section 3.2 and Figure 1 now all carry the ± 0.387.
- Feature lists are built with `sorted(list(set(...)))` throughout. The `sorted()` is load-bearing: without it, Python's set iteration order varies between interpreter sessions and silently changes XGBoost column order, which changes results. Do not "simplify" it away.

## License

MIT (see `LICENSE`). Cite the manuscript above if you use this code or reuse this validation framework.
