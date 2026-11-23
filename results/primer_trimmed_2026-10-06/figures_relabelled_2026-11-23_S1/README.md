# Figure S1 with its legend at lower left (2026-11-23)

This folder holds Figure S1 (Figure S2 in the manuscript) drawn again, with the tier names of `figures_relabelled_2026-11-23/` and its legend moved from the lower right to the lower left. At the lower right, the lower whisker of the Illumina-only tier-4 bar (−5.277 − 11.057 = −16.3) ran behind the legend. The run of record's image in `corrected_2026-11-22/figures_corrected_2026-11-22/` has the same overlap. At the lower left, no bar or whisker reaches the legend.

This drawing supersedes the Figure S1 files in `figures_relabelled_2026-11-23/`. That folder's Figure 1 is unchanged.

It was drawn in the same way as `figures_relabelled_2026-11-23/`. That folder's README gives the details:
- The script, `python/figures/figure_s1_substrate_ladder.py`, now has `loc="lower left"`.
- The inputs are the run of record's Stage 6 folder.
- It was drawn on Linux with the pinned packages, with Liberation Sans in place of Arial.

**Check.** The log lists the same input file, with the same SHA-256, and the same plotted values as the run of record's log (paths aside).

| File | What it holds |
|---|---|
| `Figure_S1_Substrate_Ladder.png`, `.pdf` | Figure S1 (manuscript Figure S2) |
| `figure_s1_substrate_ladder.log` | Everything the script printed: its input with its SHA-256, and every plotted value |
| `versions.txt` | The date, the Python and package versions, the font and the input folders |
