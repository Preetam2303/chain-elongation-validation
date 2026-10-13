# Primer-trimmed rerun: the run of record (2026-10-05 to 2026-10-13)

These files come from the laptop run of the primer-trimmed rerun (the run of record), unchanged, in the laptop's layout: the run folder `BIOTWIN_Grand_Merge\primer_trimmed_2026-10-06`, with the corrected run `corrected_2026-10-11` inside it. Every value in `results_reference.json` is read from them (or from `data/table_s2/`) by `python/build_results_reference.py`.

- Only each script's run of record is here. 03, 04, 07b, 07c, 07d and 09 were run more than once on these matrices as their code changed, and their earlier outputs stay on the laptop. So do the reads, DADA2's intermediate files and R/02's diagnostic tables.
- The files have Windows line endings, as written. The SHA-256 values in `results_reference.json` (`_about`) and in the figure logs read line endings as LF, so they match these files and a Linux run's alike.
- To rerun from these files, work on a copy. The scripts write their new dated folders inside the folders they read, and this folder should stay as the record.

| Stage | File or folder | Made by | What it holds |
|---|---|---|---|
| 1 | `stage1/dada2_overnight_log_2026-10-05.txt` | R/01b-01f | Everything the five scripts printed in their overnight run: the primer checks and each sample's read tracking. |
| 1 | `stage1/read_tracking_Duber_2020_Master.csv` | R/01a | Duber 2020's read tracking (01a was run on its own). `merged_old_run` is the same sample in the published run. |
| 2 | `ASV_Count_Matrix_Clean.csv` | R/02, R/03 | The merged, chimera-free ASV table: 121 sequencing runs x 3,672 ASVs. |
| 2 | `ASV_Taxonomy_Master.csv` | R/03 | Each ASV's SILVA v138.1 taxonomy and sequence. |
| 3 | `corrected_2026-10-11/BIOTWIN_FINAL_ML_MATRIX.csv` | R/04 | The ASV table joined to the metadata, with `data/corrections/` applied: 129 rows (117 sequencing runs), every ASV. |
| 4 | `corrected_2026-10-11/BIOTWIN_PRUNED_ML_MATRIX.csv` | 00 | The 112 modelled Illumina rows, with the 403 ASVs the 5% prevalence filter keeps. |
| 4 | `corrected_2026-10-11/BIOTWIN_GENUS_ML_MATRIX.csv` | 00 | The same rows with the ASVs summed by genus (109 genera). |
| 4 | `corrected_2026-10-11/stage4_genus_summary.csv` | 00 | Each kept genus: its ASVs, samples and share of the reads. |
| 4 | `corrected_2026-10-11/inocula_genus_profiles.csv` | 00 | The genus profiles of the 17 inoculum rows, for description only. |
| 5 | `corrected_2026-10-11/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv` | 00b | The final matrix: 112 Illumina and 54 Nanopore rows, 205 genera. |
| 5 | `corrected_2026-10-11/stage5_genus_summary.csv` | 00b | Each genus, the platforms it is found on and its share of each platform's reads. |
| 6 | `corrected_2026-10-11/stage6_corrected_2026-10-11/` | `run_stage6.py` | Run of record for 02, 05, 06, 07a, 08 and 10. |
| 6 | `corrected_2026-10-11/stage6_corrected_2026-10-11_03_04_07b_09/` | `run_stage6.py 03 04 07b 09` | Run of record for 03, 04 and 09 (its 07b output is superseded and not copied). |
| 6 | `corrected_2026-10-11/stage6_corrected_2026-10-12_07b_11/` | `run_stage6.py 07b 11` | Run of record for 07b and 11. |
| 6 | `corrected_2026-10-11/stage6_corrected_2026-10-12_07c_07d/` | `run_stage6.py 07c 07d` | Run of record for 07c and 07d. |
| 7 | `data/table_s2/` (in the repository) | 12, 12b | Table S2, with its own README. |
| 8 | `corrected_2026-10-11/figures_corrected_2026-10-13/` | `run_figures.py` | Every figure as PNG and PDF. |

Each Stage 6 folder holds each script's log (everything it printed), the tables it wrote and `versions.txt` (the date, the Python and package versions, the input files and the scripts run). The figure folder holds one log per figure script, listing each input file with its SHA-256 and every value the figure plots, and its own `versions.txt`.

Checks:
- **Against the cloud's runs.** The cloud ran the same scripts (Linux, Python 3.13.14) on its own copy of the matrices. `python python/build_results_reference.py --check` on those runs finds no value that differs from these files: 1,381 values, every one but those of Stage 1 (the `primer_trimming` section), whose files only this run has. For Table S2, the cloud ran 12b on its own run of 12 with the same BLAST file.
- **Matrices.** Cell by cell, the cloud's copies of the seven matrices hold the same values as these. They differ only in number formatting (`3.0` for `3`), in the order of tied genera in `stage4_genus_summary.csv`, and in 7 chemistry cells that carry float noise in the cloud's copies (`76.98999999999998` for `76.99`, less than 1e-9 of their size). The cloud built its R/04 matrix with a Python stand-in, because the metadata workbook R/04 reads is only on the laptop.
