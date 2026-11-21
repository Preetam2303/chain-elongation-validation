# Reactor conditions of the modelled samples

Written by `python/tables/reactor_conditions_table.py`; do not edit by hand. Rebuild with the script after any change to the matrix.

- `reactor_conditions_by_group.csv`: one row per reactor or trial group (Paper_ID + BIOREACTOR, the groups that 03b and 04 use).
- `reactor_conditions_by_group.md`: the same table, with the main columns, for reading.
- `reactor_conditions_by_sample.csv`: one row per modelled sample, with its recipe phase.
- `reactor_conditions.log`: the inputs with their SHA-256 (line endings read as LF) and the checks.

Where each column comes from:
- From the matrix: days sampled, samples, sampling units (the same-day replicate runs of Duber 2025 and Prusak 2026 form one unit; every other row is its own), `Operation_Mode`, pH, temperature, HRT and the two feed descriptors. A value that changes within a group is followed by the days it holds; in a group sampled on one day (the Duber 2022 bottles), each sample's value is listed.
- From `data/reactor_conditions/trial_recipes.csv`: the feed recipe of each phase that contains a sampled day, as the source paper gives it, with its source.
- From `data/reactor_conditions/reactor_groups.csv`: the reactor's name in the paper, the vessel and mode, whether its pH is a set point or a measured value, and the sources.

pH follows the definition of 25 Oct 2026: the culture pH at sampling, that is the set point where pH was controlled and the measured value where it was not.

`Operation_Mode` is the matrix's label. No script uses it: it is not a model input. The vessel and mode as the paper describes them are in `vessel_and_mode`.

## Log

```
Matrix: results/primer_trimmed_2026-10-06/corrected_2026-11-19/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv
sha256 98ffb7f1a5b9235fa7d51baa027b44fd922bc72bf9e99226e332d4e0352696c7  results/primer_trimmed_2026-10-06/corrected_2026-11-19/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv
sha256 24565e8d1b297d3ef90bcdcd2c16d0507f7279a96ad416822d539170c9dfe002  data/reactor_conditions/trial_recipes.csv
sha256 7e933df685dfa48e0f03b5de48e248e2e63e6c5fbe860564d10ba9d557405c40  data/reactor_conditions/reactor_groups.csv
Rows: 166; reactor or trial groups: 23; sampling units: 96
Rows by study: Brodowski 2022 8, Brodowski 2022 [Ext.] 14, Brodowski 2025 14, Duber 2020 7, Duber 2022 18, Duber 2025 51, Prusak 2026 54
Every sampled day falls in exactly one recipe phase of its group.
Recipe phases with no sample: 0
Feed descriptors that vary within a group: Brodowski 2022 [Ext.] B1 (Primary_Carbon_Signature)
Feed_Complexity: values per study: Brodowski 2022 1, Brodowski 2022 [Ext.] 1, Brodowski 2025 1, Duber 2020 1, Duber 2022 1, Duber 2025 1, Prusak 2026 1
Primary_Carbon_Signature: values per study: Brodowski 2022 1, Brodowski 2022 [Ext.] 2, Brodowski 2025 2, Duber 2020 1, Duber 2022 3, Duber 2025 1, Prusak 2026 1
```
