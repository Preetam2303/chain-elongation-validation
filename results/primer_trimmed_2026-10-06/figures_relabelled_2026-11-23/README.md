# Figures 1 and S1 with corrected labels (2026-11-23)

**Figure S1 here is superseded by `figures_relabelled_2026-11-23_S1/`**, drawn the same day with the same tier names and the legend moved to the lower left. Figure 1 here is current.

On 23 Nov 2026 Preetam asked for two figure labels to be corrected from the saved results, with no rerun. This folder holds the two figures drawn again with the corrected labels. Only the label text differs from the run of record's figures in `corrected_2026-11-22/figures_corrected_2026-11-22/`, which stay unchanged. The figures were made in the repository, not on the laptop.

- **Figure 1.** The bar "GroupKFold by Vessel (CSTRs Only)" now reads "GroupKFold by Vessel (Continuous-Reactor Subset)".
  - The label is on two lines, so the y-axis labels keep their width.
  - 04 chooses the subset by reactor ID (`BIOREACTOR` starting with B). In the run of record's final matrix it holds 145 samples from 11 reactors.
  - These include Duber 2020's reactor, which the matrix labels CSTR and the paper describes as a UASB reactor.
  - They also include the six samples of Duber 2025 R1 taken during its batch phase (days 168 and 172).
- **Figure S1** (Figure S2 in the manuscript). Each tier is named by the predictors it adds, as in Table 4. The old names were "Ops + Genera", "+ Foundational Feeds", "+ Intermediates" and "+ Downstream Proxies".
  - The names were checked against 02's tier definitions in `python/pipeline/pipeline_settings.py`.
  - "Feed descriptors" are `Feed_Complexity` and `Primary_Carbon_Signature`.

How they were drawn:
- Scripts: `python/figures/figure1_r2_comparison.py` and `python/figures/figure_s1_substrate_ladder.py`, as changed on 2026-11-23.
- Inputs: the run of record's Stage 6 folder (`corrected_2026-11-22/stage6_corrected_2026-11-22/`), set as `run_figures.py` sets it.
- Machine: Linux, with the pinned packages listed in `versions.txt`.
- Font: Arial is not installed on that machine, so the font was set to Liberation Sans after each script's import. Liberation Sans has Arial's character widths, so the layout is the same as on the laptop.

Checks:
- Each log lists the same input files, with the same SHA-256, and the same plotted values as the run of record's log (paths aside).
- The only difference is Figure 1's label column, where the log prints the line break as `\n`.

| File | What it holds |
|---|---|
| `Figure_1_Validation_Architectures.png`, `.pdf` | Figure 1 |
| `Figure_S1_Substrate_Ladder.png`, `.pdf` | Figure S1 (manuscript Figure S2) |
| `figure1_r2_comparison.log`, `figure_s1_substrate_ladder.log` | Everything each script printed: its inputs with their SHA-256, and every plotted value |
| `versions.txt` | The date, the Python and package versions, the font and the input folders |
