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
  00_prevalence_filter_genus.py      Rerun Stage 4: inocula out, spike-in out,
                                      5% prevalence filter, ASV -> genus collapse
  00b_nanopore_merge.py              Rerun Stage 5: adds the Nanopore rows
  01_genus_aggregation.py            ASV -> genus collapse (Methods 2.4; published run)
  pipeline_settings.py               Shared predictor sets, model settings and helpers
                                      for 02-10 (decision of 2026-10-10)
  02_loso_by_study.py                Leave-one-study-out: Table 4 (both columns),
                                      Table 3b's two LOSO rows (the headline),
                                      chemistry-only LOSO, training-mean benchmark
  03_permutation_test_999.py         999-shuffle label-permutation control of the
                                      headline (shuffled by sampling day)
  04_vessel_grouped_deployability.py GroupKFold by physical reactor (Table 3b)
  05_intrastudy_transferability_table5.py   Reactor-pair transfer with named pairs
                                             (Table 5; Figure 3, stage 3)
  06_cross_platform_transfer.py      Illumina -> Nanopore transfer (Table 3b)
  07a_algorithm_comparison_lightgbm_rf.py        Algorithm selection (Methods 2.7);
  07b_algorithm_comparison_lightgbm_standard.py  the target-leakage demonstration
                                                  behind excluding downstream co-products
  07c_figure3_intermediate_unscaled_xgboost.py   Figure 3, stage 2
  07d_figure3_stage1_random_forest.py            Figure 3, stage 1
  08_shared_genus_transfer_control.py  Shared-genus transfer control (Table S6)
  09_guild_membership_shap_ranks.py    SHAP rank stability across folds and the
                                        within-study Caproiciproducens check
  10_table3a_naive_splits.py           Naive, ungrouped schemes (Table 3a)
  run_stage6.py                        Runs 02-10 on a corrected run's matrices

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

**Shortcut**: `data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv` is the published run's final matrix -- if you just want to regenerate the published figures or tables, check out commit `844f478` and skip to step 4.

**Full pipeline from raw sequences:**
1. Install R packages (see `R_packages.txt`) and Python packages (`pip install -r requirements.txt`), using the **exact pinned versions**, not just compatible ones — see the reproducibility note below.
2. Run `R/` in numeric order to go from raw SRA accessions (Table 1) to the merged, pre-genus-collapse matrix. Each script has an "EDIT THIS" marker where you need to set your own local working directory -- these manage multi-stage SRA download and DADA2 output and are not one-command reruns.
   For the primer-trimmed rerun, `R/04_metadata_merge.R` writes the merged matrix, with every file in `data/corrections/` applied, into a new dated `corrected_<date>` folder inside the `primer_trimmed_*` run folder. Nothing written earlier is changed.
   `python/pipeline/00_prevalence_filter_genus.py` then leaves the 17 inoculum rows (`data/corrections/descriptive_only_samples.csv`) out of the modelling matrices and writes their genus profiles to `inocula_genus_profiles.csv` for descriptive use. It removes the ZymoBIOMICS spike-in genera (Imtechella, Allobacillus) of the Brodowski 2022 batch runs, applies the 5% prevalence filter to the 112 modelled Illumina samples (each counted once, so an ASV must be present in at least 6), and sums the kept ASVs by genus. It writes `BIOTWIN_PRUNED_ML_MATRIX.csv`, `BIOTWIN_GENUS_ML_MATRIX.csv` and `stage4_genus_summary.csv` into the same `corrected_<date>` folder.
   `python/pipeline/00b_nanopore_merge.py` then adds the 54 Nanopore rows by an outer join on genus columns. They are taken unchanged from the published matrix apart from the Hanna_2025 corrections. It writes `BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv` and `stage5_genus_summary.csv` into the same folder.
   `python/pipeline/run_stage6.py` runs scripts 02 to 10 on those matrices (Stage 6). It saves each script's printed output and the tables the scripts write to a new `stage6_corrected_<date>` folder inside `corrected_<date>`; `data/` and `results/` are not touched. The runner points each script at the corrected matrices through the environment variables `BIOTWIN_MATRIX`, `BIOTWIN_PRUNED`, `BIOTWIN_GENUS` and `BIOTWIN_RESULTS_DIR`. Run on its own, a script reads the published run's matrices in `data/` instead, with the inoculum and excluded samples left out. That gives neither the published nor the corrected numbers: the corrected numbers come from `run_stage6.py`, and the published ones from the code at commit `844f478`.
3. Run `python/pipeline/01_genus_aggregation.py` (using `data/historical/BIOTWIN_FINAL_ML_MATRIX.csv` and `ASV_Taxonomy_Master.csv`) to reproduce `data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv`, the published run's 184-sample x 170-genus matrix.
4. For the published Tables 3a, 3b, 4 and 5, run the remaining `python/pipeline/` scripts at commit `844f478` in numeric order; scripts 07a-07d use the historical snapshots in `data/historical/` (see each script's header). For the corrected run, use `run_stage6.py` (step 2).
5. Run any `python/figures/` script independently to regenerate that figure -- run from inside `python/figures/` or `python/pipeline/`, since paths are relative to each script's own location.
6. Check what you get against `results_reference.json`, which lists every reported value alongside the tolerance you should expect it to reproduce within.

`python/pipeline/03_permutation_test_999.py` refits the 7-fold LOSO 999 times and saves every shuffle's result to `permutation_null_distribution_999.csv` in the results folder. `figure2_permutation_histogram.py` still draws the published histogram from `results/`; it and the other figure scripts in `python/figures/` keep the published settings until the figures are redone for the corrected run.

## A note on exact reproducibility

Since the corrected settings of 2026-10-10, every XGBoost model in `python/pipeline/` uses the same settings (150 trees, learning rate 0.05, depth 3, `random_state=42`) with **no row or column subsampling** (`pipeline_settings.py`). With the pinned package versions, every number therefore reproduces exactly on any machine. Two models still sample rows, both seeded (`random_state=42`): 07a's LightGBM Random-Forest mode, which requires bagging, and 07d's scikit-learn random forest (Figure 3, stage 1), which draws bootstrap samples. The figure scripts in `python/figures/` have not been changed yet and still use the published settings.

The permutation test's *p*-value is a Monte-Carlo estimate over 999 shuffles. Shuffle *i* uses `random_state=i`, so the null distribution reproduces exactly.

The published numbers were made with earlier settings, in which most models used `subsample=0.8, colsample_bytree=0.8` and were therefore machine-dependent, and with a different predictor set. They reproduce with the code at commit `844f478` (see `results_reference.json` for their tolerance classes).

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
  BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv   The published run's 184-sample matrix.
                                              The default input of the scripts
                                              except 07a/07b/07c/07d below.
  study_provenance.csv                       One row per short-read study: Paper_ID,
                                              R script, paper, sequence accession used
                                              and the one the paper states, sample counts.
  corrections/
    excluded_samples.csv                     Samples removed from the data entirely,
                                              each with its reason and source.
    descriptive_only_samples.csv             Inoculum samples kept for description
                                              but never modelled (00 leaves them out).
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

`chimera_fold_genus_check.R` then compared the genus picture of those three studies at 1.5× and 8× (`top_genera_fold_1.5_vs_8.csv`, mean per-sample relative abundance). The genus picture does not change. In each study the top two genera keep their order, no named top-10 genus moves more than two places, and none changes by more than 4 percentage points. Caproiciproducens stays first or second in all three studies: 29.9% vs 31.7% (external acetate), 28.8% vs 28.9% (Brodowski 2022) and 28.2% vs 29.4% (Duber 2025). The largest shifts are Acinetobacter in the external-acetate study (36.4% to 32.7%) and sequences with no genus assigned, which gain 0.9-1.8 points at 8× (top unclassified bin in each study).

## Primer-trimmed rerun: corrected data and settings (decided 2026-10-10)

- **Modelled rows.** Every inoculum sample is left out of the models: the four shared inoculum runs and the Duber 2020 seed sludge ERR3200162, 17 rows in all. Each has caproate 0 and no chemistry measured on the inoculum itself, so a model could learn "all blank means 0". They are kept for description only. With ERR3200169 excluded as well, the models see 112 Illumina rows and 54 Nanopore rows.
- **Predictors.** Every validation scheme (Tables 3a, 3b, 4, 5 and S6) uses the same core set: pH, temperature, HRT, lactate, acetate, ethanol, butyrate, propionate, `Feed_Complexity`, `Primary_Carbon_Signature` and the genera (ASVs in the ASV-level rows).
  - Valerate, isovalerate, caprylate, heptanoate and isocaproate were measured only in the Illumina studies. They enter only Table 4's Illumina-only column and 07a's "Maximum" scope.
  - Isobutyrate, the alcohols, succinate, NaOH, lactose and DAY are never predictors in these schemes.
  - The first two stages of Figure 3 (07d, 07c) keep their broader predictor lists, DAY and the alcohols included, on purpose: they show what the earlier, looser setups give.
  - A compound that was not measured stays blank, never 0.
- **Headline.** The headline is the full-matrix leave-one-study-out result with exactly the core set. Table 4's full-matrix tier 3 and Table 3b's LOSO row are the same number. Tables 3b and 4 (02, 04, 06) also report a training-mean benchmark (predict the training rows' mean caproate), with RMSE next to R².
- **Same-day replicates.** The triplicates of Duber 2024 and the Nanopore study share one chemistry value per day. The permutation test therefore shuffles caproate by sampling day, and correlations are computed on sampling-day means.
- **Data corrections** (`data/corrections/metadata_corrections.csv`):
  - The Duber 2022 chemistry is reconverted from the author file with the formula used for every other study, and the compounds that file has are added.
  - On the Nanopore phase-boundary sampling days, pH, temperature and HRT are those of the phase that was ending.

### Genus counts: rerun vs. the earlier primer-sensitivity check

The primer-trimmed rerun kept 107 genera on 129 Illumina rows (Stage 4 run of 2026-10-10, with the published row set): 410 ASVs carrying 88.6% of the reads. The earlier sensitivity check, which stripped primer sequences with a text search instead of trimming the reads before DADA2, had about 92. The published run had 65.

- Most of the extra genera are rare ones that now cross the 5% prevalence threshold. The old run lost reads to chimeras that were really primer artefacts, so these genera fell just below the threshold.
- Two old genera are gone. Orrella and Paralcaligenes were exact two-parent chimeras of new ASVs, so they were artefacts of the old run.
- Shimwellia is a renaming: its sequence is now classified as Enterobacteriaceae without a genus.

The corrected run's own counts are in its `stage4_genus_summary.csv`.

## Resolved during development (kept here for the record)

- Figure 3's full three-stage progression (published: R²=0.217 → −0.228 → −5.869) is backed end to end: `07d_figure3_stage1_random_forest.py`, `07c_figure3_intermediate_unscaled_xgboost.py`, and `05_intrastudy_transferability_table5.py` respectively.
- The manuscript's description of the middle stage has been corrected. It previously read "XGBoost without downstream co-products", which the code contradicts: `07c` drops only sparse/junk columns and leaves Valerate, Isobutyrate and Heptanoate available as predictors — all three appear in that script's own SHAP output, alongside DAY as the single dominant feature. Section 3.5 now describes the stage accurately.
- Table 3b previously reported the headline LOSO result as R² = −0.095 with no standard deviation, while Table 4 reported the same tier-4 value as −0.095 ± 0.387. Table 3b, Section 3.2 and Figure 1 now all carry the ± 0.387.
- Feature lists are built with `sorted(list(set(...)))` throughout. The `sorted()` is load-bearing: without it, Python's set iteration order varies between interpreter sessions and silently changes XGBoost column order, which changes results. Do not "simplify" it away.

## License

MIT (see `LICENSE`). Cite the manuscript above if you use this code or reuse this validation framework.
