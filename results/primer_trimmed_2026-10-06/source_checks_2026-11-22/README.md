# Source checks of the run of record (2026-11-22)

Two checks of the final matrix against the authors' own files, made for Preetam's points 5 and 6 of 21 Nov 2026. Both were run in the repository on 2026-11-22 against the run of record, `corrected_2026-11-22/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv`. On 21 Nov the same checks against `corrected_2026-11-19/` gave the same two tables (line endings aside), because the feed-label fix of PR #32 changed no genus or caproate value. The authors' files are not in this repository; each output lists the files it read with their SHA-256.

## Nanopore genus columns (point 5)

`python python/checks/nanopore_genus_mapping.py <folder with the authors' files> --out <this folder>` wrote `nanopore_emu_to_matrix_genus.csv` and `nanopore_mapping_check_output.txt`.

The published run bridged the EMU names of Prusak 2026's Nanopore data (code label `Hanna_2025`) to SILVA genus names with a dictionary that is no longer available. The check rebuilds the bridge from the data alone: the EMU output (`Nanopore_Abundance.csv`, `Nanopore_Taxonomy.csv`: 286 taxa in 54 samples) and the sheet `SUCC_META` of `metadata_Expanding the product_Prusak.xlsx`, which maps the matrix's sample names to EMU codes.

- **120 of the 124 Nanopore genus columns** have the name of an EMU genus. Each equals the sum of that genus's EMU taxa in all 54 samples (largest difference 9.1e-11).
- **The other four** each have exactly one grouping of the remaining EMU taxa that reproduces them in all 54 samples, at the precision of the saved files:

| Matrix column | EMU taxa it holds | Largest difference over 54 samples |
|---|---|---|
| Escherichia-Shigella | EMU genus *Escherichia*: 1 taxon (unclassified *Escherichia coli* strain) | 0 |
| Clostridium_sensu_stricto_12 | EMU genus *Clostridium*: all 47 taxa | 8.7e-11 |
| Anaerostignum | EMU genus *Anaerotignum*: 2 taxa (*A. propionicum* DSM 1682, unclassified *A. lactatifermentans* strain) | 9.1e-13 |
| Unclassified | EMU "unassigned" (no lineage): 1 taxon | 0 |

- Together the 124 columns use each of the 286 EMU taxa exactly once. `nanopore_emu_to_matrix_genus.csv` lists every taxon with its matrix column and how that was established (same genus name, or inferred from the data as the only grouping that fits).
- **How uniqueness was shown** (the script's header has the method). For Escherichia-Shigella, Anaerostignum and Unclassified, the candidate taxa are linearly independent, with smallest singular values 9.95, 9.57 and 2.83 against the 1.5e-8 the proof needs, so no second grouping can fit. Clostridium_sensu_stricto_12 has 51 candidates but only 45 independent directions (its smallest singular value is zero to machine precision; the value printed depends on the numerical library). For it, the script tried every way to drop members and add non-members and found no other grouping that fits, even taking the column alone.
- On the Nanopore side, Clostridium_sensu_stricto_12 holds every EMU *Clostridium* taxon. The Illumina side (SILVA) splits Clostridium sensu stricto into 10 columns (1, 2, 3, 7, 8, 10, 11, 12, 13 and 16), and only cluster 12 shares its name with the Nanopore column.

## Caproate units (point 6)

`python python/checks/trace_caproate.py <folder with the authors' files> --out <this folder>` wrote `caproate_trace.csv` and `caproate_trace_output.txt`.

One caproate value per study, from the matrix row to the cell of the author's workbook that the audit of 6 Oct 2026 mapped it to. The matrix is in mM C. Five of the seven source columns are already in mM C. The Duber 2022 and Prusak 2026 columns are in g/L and are converted with mM C = g/L × 1000 / (M / nC), with M = 116.16 g/mol and nC = 6 for caproate. All seven agree within 0.51 mM C; the three Brodowski studies are whole numbers in the matrix. Only caproate was traced, so this says nothing about the units of the other metabolite columns.

| Study | Matrix row | Unit in the file | Value in the file | In mM C | Matrix |
|---|---|---|---|---|---|
| Duber 2020 | ERR3200168, day 107 | mM C | 278.880 | 278.880 | 278.88 |
| Duber 2022 | SRR18079352, R2, day 22 | g/L | 0.56196 | 29.027 | 29.03 |
| Duber 2025 (`Duber_2024`) | SRR28047211, B1, day 151 | mM C | 170.300 | 170.300 | 170.3 |
| Brodowski 2022 | SRR18962132, B2, day 65 | mM C | 136.467 | 136.467 | 136 |
| Brodowski 2022 [Ext.] | SRR13989643, B1, day 88 | mM C | 287.308 | 287.308 | 287 |
| Brodowski 2025 | SRR28044222, B2, day 54 | mM C | 232.800 | 232.800 | 233 |
| Prusak 2026 (`Hanna_2025`) | sample_B4, B1, day 96 | g/L | 6.515 | 336.519 | 336.519 |

Both scripts need openpyxl, which is not in `requirements.txt` because no analysis step reads an Excel file.
