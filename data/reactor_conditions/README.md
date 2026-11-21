# Reactor conditions: inputs from the source papers

The two inputs of `python/tables/reactor_conditions_table.py`, which builds the supplement's reactor-conditions table from a corrected run's final matrix. Every value in them was read from the source paper or the authors' files named in its source column.

## `trial_recipes.csv`
One row per feed phase of each reactor or trial group (39 rows, 23 groups): `Paper_ID`, `BIOREACTOR`, the phase's first and last day (`day_from`, `day_to`), its name in the paper, the feed recipe as the paper gives it, and the source.

- Taken from the writing thread's `trial_recipes_from_papers_2026-11-21.csv` (columns `recipe_from_paper` as `recipe`, and `source`). The recipes and day ranges were checked against the seven papers on 2026-11-21, and all agree.
- The check corrected the citation in 7 rows (and wrote "p. 2-3" as "pp. 2-3" in the 3 Brodowski 2022 batch rows):
  - Duber 2020 (4 rows): the periods are in Sec. 2.1 (Process configuration), PDF pp. 2-3. Sec. 2.2 is Analytical methods.
  - Duber 2022 (1 row): the paper's sections are not numbered. "Sec. 2" became "Materials and Methods".
  - Duber 2025 (2 rows): the ethanol and lactate spikes are in Sec. 3.3, PDF p. 7.
- The writing thread's columns on the matrix's feed labels are not copied. They belong to its feed-descriptor check.

## `reactor_groups.csv`
One row per reactor or trial group (23 rows):
- `name_in_paper`: e.g. BL and BLE for Brodowski 2025, R1 and R2 for Duber 2025 and Prusak 2026.
- `vessel_and_mode` and `vessel_source`.
- `ph_kind` and `ph_source`: whether the matrix's PH is a set point or a measured value. The rule is the definition of 25 Oct 2026: PH is the culture pH at sampling, so the set point where pH was controlled and the measured value where it was not. The sources and checks are in the pH provenance check of 2026-10-24 (all 112 Illumina rows), the Brodowski 2022 batch decision of 2026-11-19 (day 7; trial means of the triplicate bottles) and the Prusak 2026 conditions check of 2026-10-20 (all 54 Nanopore rows).
- `note`: how the samples relate to bottles or sequencing runs, where that matters for reading the table.
