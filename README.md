# Chain Elongation Validation Framework

Code accompanying:

> **Quantitative Fragility and Ecological Robustness in Machine-Learning-Predicted Chain Elongation: A Multi-Study, Multi-Platform Validation Framework**
> Preetam Banerjee. Water Supply and Bioeconomy Division, Poznań University of Technology. MSCA-LeAD, BIOTWIN/BioRef, WP6 Task 6.5.

## What this is

Microbiome-to-function machine learning models are typically validated within a single experimental system. This project asks what happens when a caproate-prediction model built on six previously published short-read (Illumina) chain-elongation studies plus one newly contributed long-read (Nanopore) dataset is validated properly instead: by study, by individual physical reactor, and across sequencing platforms.

The short version: naive validation reaches R² as high as 0.913 (Table 3a); every properly grouped scheme collapses that to at or below zero. The taxa the model identifies as core drivers, however, stay stable across every validation architecture, every algorithm tested, and even across sequencing platforms. That contrast — quantitative fragility alongside ecological robustness — is the paper's central finding.

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
                                      for 02-12 (decision of 2026-10-10)
  02_loso_by_study.py                Leave-one-study-out: Table 4 (both columns),
                                      Table 3b's two LOSO rows (the headline),
                                      chemistry-only LOSO, training-mean benchmark
  03_permutation_test_999.py         Label-permutation controls of the headline,
                                      999 shuffles each, by sampling day: across
                                      studies and within each study
  04_vessel_grouped_deployability.py GroupKFold by physical reactor (Table 3b)
  05_intrastudy_transferability_table5.py   Reactor-pair transfer with named pairs
                                             (Table 5; Figure 3, stage 3)
  06_cross_platform_transfer.py      Illumina -> Nanopore transfer (Table 3b)
  07a_algorithm_comparison_lightgbm_rf.py        Algorithm selection (Methods 2.7);
                                                  the target-leakage demonstration
                                                  behind excluding downstream co-products
  07b_algorithm_comparison_lightgbm_standard.py  XGBoost vs LightGBM (Methods 2.7): 80/20
                                                  split, shuffled 10-fold and LOSO; the two
                                                  naive schemes also with seeds 0-19
  07c_figure3_intermediate_unscaled_xgboost.py   Figure 3, stage 2
  07d_figure3_stage1_random_forest.py            Figure 3, stage 1
  08_shared_genus_transfer_control.py  Shared-genus transfer control (Table S6)
  09_guild_membership_shap_ranks.py    SHAP rank stability across folds and the
                                        within-study Caproiciproducens check
  10_table3a_naive_splits.py           Naive, ungrouped schemes (Table 3a)
  11_replicate_checks.py               Same-day replicates: how similar they are
                                        (Bray-Curtis), and the headline LOSO and
                                        Table 3a rerun with one row per sampling day
                                        (the naive schemes also with seeds 0-19)
  12_table_s2_shap_ranking.py          Table S2, steps 1-2: ASVs ranked by mean |SHAP|;
                                        the top 20 and the Caproiciproducens ASVs with
                                        90% of its reads written to a dated FASTA for BLAST
  12b_table_s2_from_blast.py           Table S2, step 4: the table built from
                                        BLAST's text file
  12c_genus_asvs_for_blast.py          The ASVs with 90% of one genus's reads, for
                                        BLAST, and their closest relatives (Oscillibacter)
  run_stage6.py                        Runs 02-12 on a corrected run's matrices

python/figures/        One script per figure, drawn from a corrected run's results
                       files (Stage 8 of the rerun)
  run_figures.py                     Draws every figure into a new dated folder
  figure_inputs.py                   Where the figures read their inputs and are saved
  all_samples_shap.py                The one XGBoost fit and TreeSHAP behind Figures 5 and S5
  figure1_r2_comparison.py           Master validation-collapse bar chart (Tables 3a, 3b)
  figure2_permutation_histogram.py   Null distributions of the permutation tests (03)
  figure3_duber_progression.py       The three-stage leakage-tightening figure (07d, 07c, 05)
  figure4_cohort_comparison.py       Illumina vs. Nanopore, per-sample data
  figure5_shap_genus.py              Top 10 genera by SHAP, one fit on all modelled samples
  figure6_attribution_vs_spread.py   SHAP rank across folds vs. number of studies (09)
  figure_s1_substrate_ladder.py      Table 4's tiers, both columns (02)
  figure_s3_intrastudy_transfer.py   Table 5's reactor pairs: R2 and rank correlation (05)
  figure_s4_duber_scatter.py         Table 5's Duber 2025 pair refit: measured vs. predicted
  figure_s5_combined_shap.py         Top 15 predictors by SHAP, the same fit as Figure 5

python/build_results_reference.py
                       Builds results_reference.json from the run of record,
                       and checks your own run against it (Stage 9)

results/
  primer_trimmed_2026-10-06/   The run of record of the primer-trimmed rerun,
                               unchanged from the laptop run: Stage 1 read
                               tracking, the Stage 2 ASV tables, the corrected
                               matrices, the Stage 6 tables and logs, and the
                               figures (README inside)
  chimera_check_2026-10-08/    The chimera checks (see "chimera removal" below)
  permutation_null_distribution_999.csv
                               The published run's 999 shuffles (across studies)
figures_output/         Where the published figure scripts (commit 844f478) wrote
                        PNG/PDF output; run_figures.py writes to the run folder

results_reference.json  The corrected run's values for the manuscript, read from
                        the run of record's files, with the files each comes
                        from and how exactly it reproduces. Start here.
```

Every figure in the manuscript, main and supplementary, has a corresponding script in `python/figures/`.

## Reproducing this

**Shortcut**: the corrected run of record is in `results/primer_trimmed_2026-10-06/`, with its matrices in `corrected_2026-10-11/`. To regenerate the corrected tables and figures without the sequence processing, copy that folder into a folder of your own, leaving out its `stage6_corrected_*` and `figures_corrected_*` folders. Set `GRAND_MERGE_DIR` in `run_stage6.py` and `run_figures.py` to your folder, run `run_stage6.py` (step 2) and `run_figures.py` (step 5, with `RUN_OF_RECORD` naming your new Stage 6 folder), then check your run (step 6). `data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv` is the published run's final matrix: to regenerate the published figures or tables, check out commit `844f478` and skip to step 4.

**Full pipeline from raw sequences:**
1. Install R packages (see `R_packages.txt`) and Python packages (`pip install -r requirements.txt`), using the **exact pinned versions**, not just compatible ones — see the reproducibility note below.
2. Run `R/` in numeric order to go from raw SRA accessions (Table 1) to the merged, pre-genus-collapse matrix. Each script has an "EDIT THIS" marker where you need to set your own local working directory -- these manage multi-stage SRA download and DADA2 output and are not one-command reruns.
   For the primer-trimmed rerun, `R/04_metadata_merge.R` writes the merged matrix, with every file in `data/corrections/` applied, into a new dated `corrected_<date>` folder inside the `primer_trimmed_*` run folder. Nothing written earlier is changed.
   `python/pipeline/00_prevalence_filter_genus.py` then leaves the 17 inoculum rows (`data/corrections/descriptive_only_samples.csv`) out of the modelling matrices and writes their genus profiles to `inocula_genus_profiles.csv` for descriptive use. It removes the ZymoBIOMICS spike-in genera (Imtechella, Allobacillus) of the Brodowski 2022 batch runs, applies the 5% prevalence filter to the 112 modelled Illumina samples (each counted once, so an ASV must be present in at least 6), and sums the kept ASVs by genus. It writes `BIOTWIN_PRUNED_ML_MATRIX.csv`, `BIOTWIN_GENUS_ML_MATRIX.csv` and `stage4_genus_summary.csv` into the same `corrected_<date>` folder.
   `python/pipeline/00b_nanopore_merge.py` then adds the 54 Nanopore rows by an outer join on genus columns. They are taken unchanged from the published matrix apart from the Hanna_2025 corrections. It writes `BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv` and `stage5_genus_summary.csv` into the same folder.
   `python/pipeline/run_stage6.py` runs scripts 02 to 11 on those matrices (Stage 6), and 12 (below). It saves each script's printed output and the tables the scripts write to a new `stage6_corrected_<date>` folder inside `corrected_<date>` (with the script numbers added, e.g. `_03_07b`, when only some scripts run); `data/` and `results/` are not touched. The runner points each script at the corrected matrices through the environment variables `BIOTWIN_MATRIX`, `BIOTWIN_PRUNED`, `BIOTWIN_GENUS`, `BIOTWIN_ASV_ALL` (04's `BIOTWIN_FINAL_ML_MATRIX.csv`, every ASV), `BIOTWIN_TAXONOMY` (the run folder's `ASV_Taxonomy_Master.csv`) and `BIOTWIN_RESULTS_DIR`. Run on its own, a script reads the published run's matrices in `data/` instead, with the inoculum and excluded samples left out. That gives neither the published nor the corrected numbers: the corrected numbers come from `run_stage6.py`, and the published ones from the code at commit `844f478`.
3. Run `python/pipeline/01_genus_aggregation.py` (using `data/historical/BIOTWIN_FINAL_ML_MATRIX.csv` and `ASV_Taxonomy_Master.csv`) to reproduce `data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv`, the published run's 184-sample x 170-genus matrix.
4. For the published Tables 3a, 3b, 4 and 5, run the remaining `python/pipeline/` scripts at commit `844f478` in numeric order; scripts 07a-07d use the historical snapshots in `data/historical/` (see each script's header). For the corrected run, use `run_stage6.py` (step 2).
5. `python/figures/run_figures.py` draws the figures (Stage 8) from the corrected run: the genus matrix from 00b and the tables that `run_stage6.py` wrote, never from numbers typed into a script. `RUN_OF_RECORD` at its top names the `stage6_corrected_*` folder that holds each Stage 6 script's run of record. Figure 3 reads the tables that 07c and 07d save (since 2026-10-12), so run `run_stage6.py 07c 07d` once first if your Stage 6 folders are older. `python python/figures/run_figures.py` draws every figure, and `python python/figures/run_figures.py 2 6 S4` only those named. Each figure goes as PNG and PDF into a new `figures_corrected_<date>` folder inside `corrected_<date>` (with the figure names added when only some are drawn), with one log per figure script and `versions.txt`. A log lists each input file with its SHA-256 (line endings read as LF, so Windows and Linux agree) and every value the figure plots. A figure script run on its own stops and asks for `run_figures.py`. The published figures come from the code at commit `844f478`.
6. Check what you get against `results_reference.json`. `python python/build_results_reference.py --check <your corrected_<date> folder>` reads the same values from your run (its matrices, Stage 6 folders and figure folder, and the Stage 1 and 2 files in the folder above it) and lists every value that differs from the file by more than half a unit in its last decimal; it exits with an error if any does. A section whose files your run does not have is listed as not checked. Each script's log must be in only one `stage6_corrected_*` folder, and each figure's in only one `figures_corrected_*` folder, or the check stops and says which. `python python/build_results_reference.py` rebuilds `results_reference.json` itself from the run of record in `results/primer_trimmed_2026-10-06/`.

`python/pipeline/03_permutation_test_999.py` refits the 7-fold LOSO 999 times for each of its two shuffles. It saves every shuffle's result to `permutation_null_distribution_999.csv` (across studies) and `permutation_null_within_study_999.csv` (within each study), and the *p*-values to `03_permutation_tests.csv` and, per held-out study, `03_permutation_tests_by_study.csv`, all in the results folder. `figure2_permutation_histogram.py` draws both null distributions from those files. It recomputes each count and *p*-value from every shuffle and stops if one differs from `03_permutation_tests.csv`.

**Table S2** (Stage 7 of the rerun) is rebuilt from scratch, because the published table's ASV list cannot be traced to a saved ranking. `run_stage6.py 12` ranks the ASVs by mean |SHAP| in one XGBoost fit on all 112 modelled Illumina samples (core set + ASVs as relative abundance, as Figure 5 does for genera), with each ASV's rank in the six leave-one-study-out fits as a stability check, and the number of studies each ASV is found in, as Figure 6 gives for genera. It writes `12_asv_shap_ranking.csv` (every kept ASV), `12_table_s2_asvs.csv` and `Table_S2_ASVs_<date>.fasta`. The FASTA holds the top 20 primer-free sequences and also the Caproiciproducens ASVs that together hold 90% of that genus's reads (counted over all of its ASVs in the modelled samples, before the prevalence filter), even outside the top 20, so the strain-level evidence rests on the whole genus. BLAST that FASTA by hand (Nucleotide BLAST, database "16S ribosomal RNA sequences (Bacteria and Archaea)", megablast) and save Download All -> Text. Then `python python/pipeline/12b_table_s2_from_blast.py <that folder> <BLAST text file>` writes `TableS2_<date>.csv` from the BLAST file directly: for each ASV the closest relative (every species tied with the top hit is listed), query cover, identity and accession, and for Caproiciproducens the share of the genus's reads by closest relative. A top hit of 99% or more is the closest known relative, not a species identification. The ranking, the FASTA, the BLAST file and the table from the run of record (2026-10-12, BLAST RID CTJH4F4E014) are in `data/table_s2/`, with a README listing each file and column.

`python python/pipeline/12c_genus_asvs_for_blast.py <Genus> <folder>` writes the same kind of FASTA for any other genus: the ASVs that together hold 90% of its reads in the modelled Illumina samples, read from the run of record in `results/primer_trimmed_2026-10-06/`, with each ASV's SHAP rank from 12's ranking. No model is fitted, so it gives the same file on every machine. Run again with the BLAST text file as a third argument, it adds the closest relatives as 12b does. `data/oscillibacter_blast/` holds it for Oscillibacter (2026-10-14), whose identity the Discussion needs from the corrected run's primer-free sequences.

## A note on exact reproducibility

Since the corrected settings of 2026-10-10, every XGBoost model in `python/pipeline/` uses the same settings (150 trees, learning rate 0.05, depth 3, `random_state=42`) with **no row or column subsampling** (`pipeline_settings.py`). With the pinned package versions, every number therefore reproduces exactly on any machine. The cloud's own runs of the same scripts (Linux, Python 3.13.14, from its own copy of the matrices) agree with the laptop run of record (Windows 11, Python 3.13.9) on all 1,523 values of `results_reference.json` that both runs have, which is every value but those of Stage 1 (`--check`, step 6). Two models still sample rows, both seeded (`random_state=42`): 07a's LightGBM Random-Forest mode, which requires bagging, and 07d's scikit-learn random forest (Figure 3, stage 1), which draws bootstrap samples. The figure scripts in `python/figures/` use the same settings. Figures 5 and S5 come from one XGBoost fit on all 166 modelled samples with TreeSHAP on the same samples (no split, so no seed). Figure S4 refits Table 5's Duber 2025 pair exactly as 05 does and stops if its R², RMSE or rank correlation differs from 05's table. The only random parts of the figures are visual (the jitter of the points in Figures 4, 5 and 6, and the order in which Figure S5 draws its dots), and they are seeded.

The permutation tests' *p*-values are Monte-Carlo estimates over 999 shuffles each. Shuffle *i* uses `random_state=i` (across studies) or `np.random.RandomState(i)` (within each study), so both null distributions reproduce exactly.

The vessel folds of 04 and 09 are scikit-learn's GroupKFold folds, with ties between vessels of the same size broken by name (`group_kfold_splits` in `pipeline_settings.py`). GroupKFold itself orders such ties with numpy's default sort, whose order differs between processors. Before 2026-10-11, the laptop and the cloud therefore put some of these vessels, such as the two Nanopore reactors (27 rows each), into different folds and got different vessel-grouped numbers.

The published numbers were made with earlier settings, in which most models used `subsample=0.8, colsample_bytree=0.8` and were therefore machine-dependent, and with a different predictor set. They reproduce with the code at commit `844f478`. Their values and tolerance classes are in `results_reference.json` as it stood at tag `draft-v1-2026-10-04` (commit `29635f7`).

**Environment used for every reported number (the laptop run of record):**

```
Python 3.13.9
Windows 11
package versions exactly as pinned in requirements.txt and R_packages.txt
```

## Data availability

Raw sequencing data for the six short-read studies are available under their original NCBI/ENA accessions (Table 1 of the manuscript). The Nanopore dataset (Prusak et al., 2026) is available under BioProject PRJNA1222907.

## Data

```
data/
  BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv   The published run's 184-sample matrix.
                                              The default input of the scripts
                                              except 07a/07c/07d below (07b reads
                                              it and the pruned matrix).
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
                                              genus aggregation; 07b's 80/20 split),
                                              and by 10 and 11 (Table 3a's
                                              ASV-level schemes).
    BIOTWIN_FINAL_ML_MATRIX.csv,
    ASV_Taxonomy_Master.csv                  Inputs to 01_genus_aggregation.py.
  table_s2/                                  Table S2 of the rerun (2026-10-12): 12's
                                              ASV ranking and FASTA, the BLAST text
                                              file and 12b's TableS2_2026-10-12.csv,
                                              unchanged from the laptop run (README).
  oscillibacter_blast/                       12c's Oscillibacter ASVs and FASTA (2026-10-14).
```

All four dataset snapshots are provided so the complete history behind Figure 3 -- not just its final, headline number -- stays runnable, in keeping with this repository's whole premise: every reported number traceable to a script that produces it.

## Primer-trimmed rerun: primer removal

The published run kept the 341F/805R primers in the reads. Left on the reads, each combination of the primers' wobble bases looked like a different sequence to DADA2. The rerun therefore removes both primers before DADA2: R/01a-01f cut them off every read pair with cutadapt 5.2.

```
cutadapt -g ^CCTACGGGNGGCWGCAG -G ^GACTACHVGGGTATCTAATCC -e 2 --discard-untrimmed --pair-filter=any
```

- `^` anchors each primer at the start of its read, and `-e 2` allows up to 2 mismatches or indels. A fixed `trimLeft` cut is not used, because some primer copies are 1-2 bases shorter than 17 or 21 bases; cutadapt finds where the primer ends in each read.
- A read pair in which either primer is not found is discarded. The amplicon (about 464 bp) is longer than the reads, so the opposite primer never appears at the read ends.
- `truncLen=c(233,229)` is the old 250 minus the primer lengths (17, 21), so the reads cover the same stretch of the 16S gene as before.
- The reverse primer is the one some papers call 785R.

Read pairs in which cutadapt found both primers (the only pairs DADA2 sees), from each script's read tracking in `results/primer_trimmed_2026-10-06/stage1/` (`results_reference.json`, `primer_trimming`):

| Study | Samples | Read pairs in | With both primers | Lowest-highest sample |
|---|---|---|---|---|
| Duber 2020 | 9 | 938,025 | 95.9% | 94.6-96.4% |
| Duber 2022 | 19 | 1,831,322 | 97.6% | 97.2-98.2% |
| Brodowski 2022 | 10 | 1,281,663 | 98.4% | 98.0-98.9% |
| Brodowski 2022 [Ext.] | 16 | 1,957,049 | 98.4% | 98.1-98.7% |
| Brodowski 2025 | 15 | 2,297,696 | 97.5% | 97.3-97.6% |
| Duber 2025 (`Duber_2024`) | 53 | 4,339,860 | 97.1% | 96.4-97.7% |

Brodowski 2022 [Ext.]'s 16 samples include the mock community, which R/01d removes after DADA2. R/02 then merges the six studies and removes chimeras, and R/03 assigns taxonomy: 121 sequencing runs and 3,672 ASVs. No ASV of that table starts with 341F or ends with the reverse complement of 805R.

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
- **Same-day replicates.** The triplicates of Duber 2024 and the Nanopore study share one chemistry value per day. The permutation tests therefore shuffle caproate by sampling day, and correlations are computed on sampling-day means. Every sample stays a row in the main analysis; 11 shows how similar a day's three samples are, and reruns the headline LOSO and Table 3a's naive schemes with one row per sampling day for comparison (decided 2026-10-07).
- **Permutation tests** (decided 2026-10-11). Caproate is shuffled in two ways, reported side by side. Shuffling across studies, as published, mixes values between studies and so erases the between-study differences that make the real task hard. Shuffling within each study keeps those differences, so it asks whether the predictors carry information that transfers beyond study-level patterns. Both report R² and RMSE.
- **Algorithm comparison** (decided 2026-10-11). The ranking of XGBoost and LightGBM depends on the validation scheme, so 07b compares them on the single 80/20 split, shuffled 10-fold and LOSO.
- **Seeds of the naive schemes.** A naive score depends on which rows land in the test set, so one seed can mislead. 07b and 11 also run their naive schemes with seeds 0-19 and report the median and range; the seed-42 numbers stay as they are, so Table 3a and 07b's first table do not change.
- **Data corrections** (`data/corrections/metadata_corrections.csv`):
  - The Duber 2022 chemistry is reconverted from the author file with the formula used for every other study, and the compounds that file has are added.
  - On the Nanopore phase-boundary sampling days, pH, temperature and HRT are those of the phase that was ending.

### Genus counts: rerun vs. the earlier primer-sensitivity check

The primer-trimmed rerun kept 107 genera on 129 Illumina rows (Stage 4 run of 2026-10-10, with the published row set): 410 ASVs carrying 88.6% of the reads. The earlier sensitivity check, which stripped primer sequences with a text search instead of trimming the reads before DADA2, had about 92. The published run had 65.

- Most of the extra genera are rare ones that now cross the 5% prevalence threshold. The old run lost reads to chimeras that were really primer artefacts, so these genera fell just below the threshold.
- Two old genera are gone. Orrella and Paralcaligenes were exact two-parent chimeras of new ASVs, so they were artefacts of the old run.
- Shimwellia is a renaming: its sequence is now classified as Enterobacteriaceae without a genus.

The corrected run (2026-10-11), with the inoculum rows left out, keeps 403 ASVs carrying 90.3% of the reads of the 112 modelled Illumina samples, in 109 genera (`results/primer_trimmed_2026-10-06/corrected_2026-10-11/stage4_genus_summary.csv`; `results_reference.json`, `data_construction`).

## Resolved during development (kept here for the record)

- Figure 3's full three-stage progression (published: R²=0.217 → −0.228 → −5.869; corrected run: 0.155 → −0.709 → −5.296) is backed end to end: `07d_figure3_stage1_random_forest.py`, `07c_figure3_intermediate_unscaled_xgboost.py`, and `05_intrastudy_transferability_table5.py` respectively.
- The manuscript's description of the middle stage has been corrected. It previously read "XGBoost without downstream co-products", which the code contradicts: `07c` drops only sparse/junk columns and leaves Valerate, Isobutyrate and Heptanoate available as predictors — in the published run all three appear in that script's own SHAP output, alongside DAY as the single dominant feature. Section 3.5 now describes the stage accurately. In the corrected run, DAY still dominates, followed by Clostridium sensu stricto 11, pH and isobutyrate; valerate and heptanoate are not in the top 10 (`results_reference.json`, `figure_3_leakage_tightening_progression`).
- Table 3b previously reported the headline LOSO result as R² = −0.095 with no standard deviation, while Table 4 reported the same tier-4 value as −0.095 ± 0.387. Table 3b, Section 3.2 and Figure 1 now all carry the ± 0.387.
- Feature lists are built with `sorted(list(set(...)))` throughout. The `sorted()` is load-bearing: without it, Python's set iteration order varies between interpreter sessions and silently changes XGBoost column order, which changes results. Do not "simplify" it away.

## License

MIT (see `LICENSE`). Cite the manuscript above if you use this code or reuse this validation framework.
