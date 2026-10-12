# Table S2: ranking, FASTA, BLAST file and table (rerun, 2026-10-12)

These files come from one laptop run (the run of record) of the primer-trimmed rerun, unchanged:
- `python/pipeline/run_stage6.py 12`, written to `corrected_2026-10-11\stage6_corrected_2026-10-12_12`;
- an NCBI web BLAST of the FASTA;
- `python/pipeline/12b_table_s2_from_blast.py` on that folder and the BLAST file.

The CSV files and the log have Windows line endings. The method was decided by Preetam on 2026-10-12. See the README's **Table S2** paragraph and the two scripts' headers.

| File | Made by | What it holds |
|---|---|---|
| `12_asv_shap_ranking.csv` | 12 | Every kept ASV (403), ranked by mean \|SHAP\| in one XGBoost fit on all 112 modelled Illumina samples. Also each ASV's samples, studies, share of its genus's reads, and rank in the six leave-one-study-out fits. |
| `12_table_s2_asvs.csv` | 12 | The 36 ASVs in the FASTA, with their sequences: the top 20, plus the other Caproiciproducens ASVs that together hold 90% of the genus's reads. |
| `Table_S2_ASVs_2026-10-12.fasta` | 12 | The primer-free sequences that were BLASTed. |
| `12_table_s2_shap_ranking.log`, `versions.txt` | 12 / runner | Everything 12 printed, and the package versions and input files. |
| `CTJH4F4E014-Alignment.txt` | NCBI BLAST | RID CTJH4F4E014, saved with Download All → Text. Nucleotide BLAST (BLASTN, megablast) against "16S ribosomal RNA sequences (Bacteria and Archaea)" (rRNA_typestrains/16S_ribosomal_RNA). |
| `TableS2_2026-10-12.csv` | 12b | Table S2: one row per ASV, in the FASTA's order. |

Columns of `TableS2_2026-10-12.csv`:
- `rank`, `mean_abs_shap`: the ASV's rank and mean |SHAP| in the all-samples fit. Tied ASVs get their average rank, and an ASV with no attribution has mean |SHAP| 0.
- `included_as`: in the SHAP top 20, in the Caproiciproducens 90% set, or in both.
- `silva_taxon`, `genus`: the SILVA genus, or the lowest rank SILVA assigned.
- `n_samples`, `n_studies`: the modelled samples and studies with a count above 0 (studies as in Figure 6).
- `pct_of_genus_reads`, `cum_pct_of_genus_reads`: the ASV's share of its genus's reads over all of that genus's ASVs before the prevalence filter. The cumulative share is filled in for the Caproiciproducens 90% set only.
- `closest_relatives`, `accessions`: the top hit and every other species with the same max score, identity and query cover, all listed. A top hit of 99% or more is the closest known relative, not a species identification.
- `query_cover_pct`, `identity_pct`, `at_least_99`: the top hit's cover and identity, and whether the identity is at least 99%.
- `n_species_tied`, `next_species`, `next_identity_pct`: how many species tie for the top hit, and the best species outside the tie.
- `n_hits`, `RID`: the number of hits BLAST returned, and the BLAST search they come from.

Checks:
- **Cloud rerun of 12.** On the same matrices, it gives the same ranking, table and FASTA.
- **Cloud rerun of 12b.** On these files, it gives a byte-identical `TableS2_2026-10-12.csv`.
- **Independent read.** A second reader that uses only the BLAST file's alignment sections agrees with the table for all 36 ASVs.
