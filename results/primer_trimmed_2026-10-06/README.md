# Primer-trimmed rerun: the run of record

These files come from the laptop runs of the primer-trimmed rerun, unchanged, in the laptop's layout: the run folder `BIOTWIN_Grand_Merge\primer_trimmed_2026-10-06`, with the corrected runs inside it.

- **Run of record since 2026-11-19: `corrected_2026-11-19/`.** This is the full laptop rerun from R/04 onward after PR #30, which set the three Brodowski 2022 batch rows (SRR18962128-30) to the day-7 interpretation: PH 6.50 / 6.56 / 6.44 and DAY 7. Nothing else changed in the code or data. Every value in `results_reference.json` is read from it by `python/build_results_reference.py`. Table S2 (12 and 12b) is in its Stage 6 folder.
- **Earlier run of record (2026-10-11 to 2026-11-19): `corrected_2026-10-11/`**, kept unchanged. Its values are in `results_reference.json` as of commit 95420b1. `differences_2026-10-11_to_2026-11-19.csv` lists every value of that file next to the same value in the new run, and whether it differs (473 of 1,838 do).
- Only each script's run of record is here. The reads, DADA2's intermediate files and R/02's diagnostic tables stay on the laptop.
- The files have Windows line endings, as written. The SHA-256 values in `results_reference.json` (`_about`) and in the figure logs read line endings as LF, so they match these files and a Linux run's alike.
- To rerun from these files, work on a copy. The scripts write their new dated folders inside the folders they read, and this folder should stay as the record.

## Reactor-conditions table (supplement)

`reactor_conditions_2026-11-21/` holds the reactor-conditions table built from the run of record's final matrix (`corrected_2026-11-19/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv`) by `python/tables/reactor_conditions_table.py`: one row per reactor or trial group (23 groups, 166 samples, 96 sampling units), one row per sample, and a log with the inputs' SHA-256. It is made in the repository, not on the laptop, and changes nothing in the run folders. Its README says where each column comes from.

## Shared by both runs (Stages 1 and 2)

| Stage | File or folder | Made by | What it holds |
|---|---|---|---|
| 1 | `stage1/dada2_overnight_log_2026-10-05.txt` | R/01b-01f | Everything the five scripts printed in their overnight run: the primer checks and each sample's read tracking. |
| 1 | `stage1/read_tracking_Duber_2020_Master.csv` | R/01a | Duber 2020's read tracking (01a was run on its own). `merged_old_run` is the same sample in the published run. |
| 2 | `ASV_Count_Matrix_Clean.csv` | R/02, R/03 | The merged, chimera-free ASV table: 121 sequencing runs x 3,672 ASVs. |
| 2 | `ASV_Taxonomy_Master.csv` | R/03 | Each ASV's SILVA v138.1 taxonomy and sequence. |

## Run of record: `corrected_2026-11-19/` (laptop, 2026-11-19)

| Stage | File or folder | Made by | What it holds |
|---|---|---|---|
| 3 | `corrected_2026-11-19/BIOTWIN_FINAL_ML_MATRIX.csv` | R/04 | The ASV table joined to the metadata, with `data/corrections/` applied (339 value corrections): 129 rows (117 sequencing runs), every ASV. |
| 4 | `corrected_2026-11-19/BIOTWIN_PRUNED_ML_MATRIX.csv` | 00 | The 112 modelled Illumina rows, with the 403 ASVs the 5% prevalence filter keeps. |
| 4 | `corrected_2026-11-19/BIOTWIN_GENUS_ML_MATRIX.csv` | 00 | The same rows with the ASVs summed by genus (109 genera). |
| 4 | `corrected_2026-11-19/stage4_genus_summary.csv` | 00 | Each kept genus: its ASVs, samples and share of the reads. |
| 4 | `corrected_2026-11-19/inocula_genus_profiles.csv` | 00 | The genus profiles of the 17 inoculum rows, for description only. |
| 5 | `corrected_2026-11-19/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv` | 00b | The final matrix: 112 Illumina and 54 Nanopore rows, 205 genera. |
| 5 | `corrected_2026-11-19/stage5_genus_summary.csv` | 00b | Each genus, the platforms it is found on and its share of each platform's reads. |
| 6, 7 | `corrected_2026-11-19/stage6_corrected_2026-11-19/` | `run_stage6.py`, then 12b | Every Stage 6 script (02 to 11 and 03b) and Table S2: 12's ranking and FASTA, and 12b's `TableS2_2026-11-19.csv`, built from the BLAST file `data/table_s2/CTJH4F4E014-Alignment.txt`. 12 picked the same ASVs as on 2026-10-12, with identical sequences, so that BLAST file still applies. |
| 8 | `corrected_2026-11-19/figures_corrected_2026-11-19/` | `run_figures.py` | Every figure as PNG and PDF. The runner printed `stage6_corrected_2026-11-19` as the folder of all nine Stage 6 scripts it reads. |

## Earlier run of record: `corrected_2026-10-11/` (laptop, 2026-10-11 to 2026-10-20)

Same Stages 3 to 5 files as above, from the matrices before PR #30. 03, 04, 07b, 07c, 07d and 09 were run more than once on these matrices as their code changed; only each script's run of record is here.

| Stage | File or folder | Made by | What it holds |
|---|---|---|---|
| 6 | `corrected_2026-10-11/stage6_corrected_2026-10-11/` | `run_stage6.py` | Run of record for 02, 05, 06, 07a, 08 and 10. |
| 6 | `corrected_2026-10-11/stage6_corrected_2026-10-11_03_04_07b_09/` | `run_stage6.py 03 04 07b 09` | Run of record for 03, 04 and 09 (its 07b output is superseded and not copied). |
| 6 | `corrected_2026-10-11/stage6_corrected_2026-10-12_07b_11/` | `run_stage6.py 07b 11` | Run of record for 07b and 11. |
| 6 | `corrected_2026-10-11/stage6_corrected_2026-10-12_07c_07d/` | `run_stage6.py 07c 07d` | Run of record for 07c and 07d. |
| 6 | `corrected_2026-10-11/stage6_corrected_2026-10-20_03b/` | `run_stage6.py 03b` | Run of record for 03b (added 2026-10-19; exploratory sensitivity analyses for 03). |
| 7 | `data/table_s2/` (in the repository; 2026-10-12) | 12, 12b | Table S2, with its own README. |
| 8 | `corrected_2026-10-11/figures_corrected_2026-10-14/` | `run_figures.py` | Every figure as PNG and PDF. |

Each Stage 6 folder holds each script's log (everything it printed), the tables it wrote and `versions.txt` (the date, the Python and package versions, the input files and the scripts run). The figure folder holds one log per figure script, listing each input file with its SHA-256 and every value the figure plots, and its own `versions.txt`.

Checks of the run of record (`corrected_2026-11-19/`):
- **Against the cloud's run of the same values.** The cloud ran the same pipeline on 2026-10-25 (Linux, Python 3.13.14), with the three batch rows set to the same day-7 values. `compare_runs.py` reads the same 1,838 values from both runs, and none differs. Of the 33 tables and FASTA files, 31 are byte-identical (line endings aside). Of the other two, `BIOTWIN_FINAL_ML_MATRIX.csv` holds the same values, with one cell written `6.5` where the cloud's copy has `6.50`, and `stage4_genus_summary.csv` has two pairs of tied genera in the other order.
- **Against itself.** `python python/build_results_reference.py --check results/primer_trimmed_2026-10-06/corrected_2026-11-19` compares 1,838 values, and none differs.

Checks of the earlier run of record (`corrected_2026-10-11/`):
- **Against the cloud's runs.** The cloud ran the same scripts (Linux, Python 3.13.14) on its own copy of the matrices. `python python/build_results_reference.py --check` on those runs finds no value that differs from these files: 1,786 values, every one but those of Stage 1 (the `primer_trimming` section), whose files only this run has. 03b's four tables are also byte-identical to the cloud's (line endings aside). For Table S2, the cloud ran 12b on its own run of 12 with the same BLAST file.
- **Matrices.** Cell by cell, the cloud's copies of the seven matrices hold the same values as these. They differ only in number formatting (`3.0` for `3`), in the order of tied genera in `stage4_genus_summary.csv`, and in 7 chemistry cells that carry float noise in the cloud's copies (`76.98999999999998` for `76.99`, less than 1e-9 of their size). The cloud built its R/04 matrix with a Python stand-in, because the metadata workbook R/04 reads is only on the laptop.
