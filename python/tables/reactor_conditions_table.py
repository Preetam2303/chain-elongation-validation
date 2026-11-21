# reactor_conditions_table.py
# The reactor-conditions table of the supplement: one row per reactor or trial
# group (Paper_ID + BIOREACTOR, the 23 groups that 03b and 04 use), covering
# every modelled sample. Built by script so it can be rebuilt after any change
# to the matrix.
#
# What comes from where:
#   - the matrix (BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv of a corrected run):
#     days sampled, samples, sampling units, Operation_Mode, PH, TEMP, HRT and
#     the two feed descriptors (Feed_Complexity, Primary_Carbon_Signature);
#   - data/reactor_conditions/trial_recipes.csv: each trial's feed recipe by
#     phase, as the source paper gives it, with its source;
#   - data/reactor_conditions/reactor_groups.csv: the reactor's name in the
#     paper, the vessel and mode, whether its pH is a set point or a measured
#     value (PH = culture pH at sampling: the set point if pH was controlled,
#     the measured value if not), and the sources.
# Nothing in the table is typed in: every value is read from these three files.
#
# Usage, from the repository root:
#   python python/tables/reactor_conditions_table.py
#       reads the run of record's matrix (RUN_OF_RECORD in build_results_reference.py)
#       and writes results/primer_trimmed_2026-10-06/reactor_conditions_<today>/
#   python python/tables/reactor_conditions_table.py --matrix-dir <corrected_<date> folder> --out <new folder>
# The output folder must not exist yet; nothing is overwritten.

import argparse
import datetime
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "python"))
from build_results_reference import RUN_OF_RECORD  # noqa: E402
import pipeline_settings as ps  # noqa: E402  (on the path via build_results_reference)
from figure_inputs import STUDY_LABELS, sha256_lf  # noqa: E402

INPUTS = REPO / "data" / "reactor_conditions"
MATRIX_NAME = "BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"
EXPECTED = dict(rows=166, groups=23)


README = """# Reactor conditions of the modelled samples

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
{log}
```
"""


def day_number(day):
    """DAY as a number: the Nanopore days are written d0, d15, ..."""
    return int(str(day).strip().lstrip("d"))


def fmt(x, col):
    """pH with one or two decimals (6.0, 6.56), temperature and HRT as short as exact (30, 1.25)."""
    if col == "PH":
        text = f"{x:.2f}".rstrip("0")
        return text + "0" if text.endswith(".") else text
    return f"{x:g}"


def join_days(days):
    return ", ".join(str(d) for d in sorted(set(days)))


def days_text(days):
    return ("day " if len(set(days)) == 1 else "days ") + join_days(days)


def by_day(sub, col, blank):
    """One value; or, on a single sampling day, each sample's value (one per bottle);
    or each value followed by the days it holds, in order of first day."""
    vals = sub[col].astype(object).map(lambda v: blank if pd.isna(v) else (fmt(v, col) if isinstance(v, (int, float)) else str(v)))
    if vals.nunique() == 1:
        return vals.iloc[0]
    if sub["day"].nunique() == 1:
        return " / ".join(vals[sub[col].sort_values().index])
    first = sub.groupby(vals)["day"].min().sort_values()
    return "; ".join(f"{v} ({days_text(sub.loc[vals == v, 'day'])})" for v in first.index)


def recipe_phase(recipes, paper, reactor, day):
    hit = recipes[(recipes.Paper_ID == paper) & (recipes.BIOREACTOR == reactor)
                  & (recipes.day_from <= day) & (recipes.day_to >= day)]
    if len(hit) != 1:
        raise SystemExit(f"{paper} {reactor} day {day}: {len(hit)} recipe phases contain this day (need exactly 1).")
    return hit.index[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix-dir", type=Path, default=RUN_OF_RECORD,
                    help="corrected_<date> folder holding " + MATRIX_NAME)
    ap.add_argument("--out", type=Path, default=None,
                    help="new output folder (default results/primer_trimmed_2026-10-06/reactor_conditions_<today>)")
    a = ap.parse_args()
    out = a.out or REPO / "results" / "primer_trimmed_2026-10-06" / f"reactor_conditions_{datetime.date.today():%Y-%m-%d}"
    if out.exists():
        raise SystemExit(f"{out} already exists. Nothing is overwritten; give a new folder with --out.")

    matrix = a.matrix_dir / MATRIX_NAME
    recipes_path, groups_path = INPUTS / "trial_recipes.csv", INPUTS / "reactor_groups.csv"
    inputs = [matrix, recipes_path, groups_path]
    for p in inputs:
        print(f"Input: {p}\n       sha256 {sha256_lf(p)}")

    df = ps.load_matrix(matrix, relative=False)[0].copy()
    df["day"] = df["DAY"].map(day_number)
    df["unit"] = ps.sampling_units(df)
    df["group"] = df["Paper_ID"] + "|" + df["BIOREACTOR"]
    recipes = pd.read_csv(recipes_path)
    groups = pd.read_csv(groups_path).fillna("")
    groups["group"] = groups["Paper_ID"] + "|" + groups["BIOREACTOR"]

    # Every group in the matrix has exactly one entry in each input file, and no input names a group the matrix lacks.
    in_matrix = set(df["group"])
    for name, keys in (("reactor_groups.csv", groups["group"]),
                       ("trial_recipes.csv", recipes["Paper_ID"] + "|" + recipes["BIOREACTOR"])):
        if set(keys) != in_matrix:
            raise SystemExit(f"{name} and the matrix name different groups: "
                             f"only in matrix {sorted(in_matrix - set(keys))}, only in {name} {sorted(set(keys) - in_matrix)}")
    if groups["group"].duplicated().any():
        raise SystemExit("reactor_groups.csv lists a group twice.")
    if len(df) != EXPECTED["rows"] or len(in_matrix) != EXPECTED["groups"]:
        print(f"NOTE: {len(df)} rows in {len(in_matrix)} groups (the run of record has "
              f"{EXPECTED['rows']} rows in {EXPECTED['groups']} groups).")

    df["phase_row"] = [recipe_phase(recipes, p, r, d) for p, r, d in zip(df.Paper_ID, df.BIOREACTOR, df.day)]
    unsampled = sorted(set(recipes.index) - set(df["phase_row"]))

    # Per sample
    g = groups.set_index("group")
    rec = recipes.loc[df["phase_row"]].reset_index(drop=True)
    per_sample = pd.DataFrame({
        "study": df["Paper_ID"].map(STUDY_LABELS),
        "Paper_ID": df["Paper_ID"],
        "reactor_or_trial": df["BIOREACTOR"],
        "name_in_paper": df["group"].map(g["name_in_paper"]),
        "Sample_ID": df["Sample_ID"],
        "day": df["day"],
        "sampling_unit": df["unit"].where(df["Paper_ID"].isin(ps.REPLICATE_STUDIES), df["Sample_ID"]),
        "Operation_Mode": df["Operation_Mode"],
        "PH": df["PH"],
        "pH_kind": df["group"].map(g["ph_kind"]),
        "TEMP_C": df["TEMP"],
        "HRT_d": df["HRT"],
        "Feed_Complexity": df["Feed_Complexity"].astype(str),
        "Primary_Carbon_Signature": df["Primary_Carbon_Signature"].astype(str),
        "recipe_phase": rec["phase"],
        "recipe": rec["recipe"],
        "recipe_source": rec["source"],
    })
    order = sorted(df["group"].unique(), key=lambda k: (STUDY_LABELS[k.split("|")[0]], k.split("|")[1][0],
                                                         int(k.split("|")[1][1:])))
    per_sample["_g"] = pd.Categorical(df["group"], order, ordered=True)
    per_sample = per_sample.sort_values(["_g", "day", "Sample_ID"]).drop(columns="_g")

    # Per group
    rows = []
    for k in order:
        sub = df[df["group"] == k].sort_values("day")
        info = g.loc[k]
        phases = []
        for i in sorted(sub["phase_row"].unique(), key=lambda i: recipes.loc[i, "day_from"]):
            r = recipes.loc[i]
            phases.append(f"{r.phase} (days {r.day_from}-{r.day_to}; sampled {join_days(sub.loc[sub.phase_row == i, 'day'])}): {r.recipe}")
        sources = list(dict.fromkeys(recipes.loc[sorted(sub["phase_row"].unique()), "source"]))
        rows.append({
            "study": STUDY_LABELS[sub["Paper_ID"].iloc[0]],
            "Paper_ID": sub["Paper_ID"].iloc[0],
            "reactor_or_trial": sub["BIOREACTOR"].iloc[0],
            "name_in_paper": info["name_in_paper"],
            "vessel_and_mode": info["vessel_and_mode"],
            "vessel_source": info["vessel_source"],
            "Operation_Mode_in_matrix": by_day(sub, "Operation_Mode", ""),
            "days_sampled": join_days(sub["day"]),
            "samples": len(sub),
            "sampling_units": sub["unit"].nunique(),
            "feed_recipe_by_phase": " | ".join(phases),
            "recipe_source": " | ".join(sources),
            "Feed_Complexity_in_matrix": by_day(sub, "Feed_Complexity", ""),
            "Primary_Carbon_Signature_in_matrix": by_day(sub, "Primary_Carbon_Signature", ""),
            "pH": by_day(sub, "PH", "blank"),
            "pH_kind": info["ph_kind"],
            "pH_source": info["ph_source"],
            "temperature_C": by_day(sub, "TEMP", "blank"),
            "HRT_d": by_day(sub, "HRT", "none (batch)"),
            "note": info["note"],
        })
    per_group = pd.DataFrame(rows)

    # Checks
    assert per_group["samples"].sum() == len(df) == len(per_sample)
    assert per_group["sampling_units"].sum() == df["unit"].nunique()
    assert len(per_group) == len(in_matrix)

    out.mkdir(parents=True)
    per_group.to_csv(out / "reactor_conditions_by_group.csv", index=False, lineterminator="\n")
    per_sample.to_csv(out / "reactor_conditions_by_sample.csv", index=False, lineterminator="\n")

    md = ["| Study | Reactor or trial | Vessel and mode | Days sampled | Samples (sampling units) | Feed recipe by phase | "
          "Feed descriptors in the matrix | pH | pH kind | Temperature (°C) | HRT (d) |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in per_group.itertuples():
        name = r.reactor_or_trial if r.name_in_paper == r.reactor_or_trial else f"{r.reactor_or_trial} ({r.name_in_paper})"
        md.append(f"| {r.study} | {name} | {r.vessel_and_mode} | {r.days_sampled} | {r.samples} ({r.sampling_units}) | "
                  f"{r.feed_recipe_by_phase.replace(' | ', '<br>')} | {r.Feed_Complexity_in_matrix}; {r.Primary_Carbon_Signature_in_matrix} | "
                  f"{r.pH} | {r.pH_kind} | {r.temperature_C} | {r.HRT_d} |")
    (out / "reactor_conditions_by_group.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    log = [f"Matrix: {matrix.relative_to(REPO) if matrix.is_relative_to(REPO) else matrix}",
           *[f"sha256 {sha256_lf(p)}  {p.relative_to(REPO) if p.is_relative_to(REPO) else p}" for p in inputs],
           f"Rows: {len(df)}; reactor or trial groups: {len(per_group)}; sampling units: {df['unit'].nunique()}",
           f"Rows by study: " + ", ".join(f"{STUDY_LABELS[s]} {n}" for s, n in df.groupby('Paper_ID').size().items()),
           "Every sampled day falls in exactly one recipe phase of its group.",
           f"Recipe phases with no sample: {len(unsampled)}" + (": " + "; ".join(
               f"{recipes.loc[i, 'Paper_ID']} {recipes.loc[i, 'BIOREACTOR']} {recipes.loc[i, 'phase']}" for i in unsampled)
               if unsampled else ""),
           "Feed descriptors that vary within a group: " + (", ".join(
               f"{r.study} {r.reactor_or_trial} ({c})" for r in per_group.itertuples()
               for c, v in (("Feed_Complexity", r.Feed_Complexity_in_matrix),
                            ("Primary_Carbon_Signature", r.Primary_Carbon_Signature_in_matrix)) if "(day" in v) or "none"),
           *[f"{c}: values per study: " + ", ".join(f"{STUDY_LABELS[s]} {n}" for s, n in df.groupby("Paper_ID")[c].nunique().items())
             for c in ps.CAT_COLS]]
    (out / "reactor_conditions.log").write_text("\n".join(log) + "\n", encoding="utf-8")
    (out / "README.md").write_text(README.format(log="\n".join(log)), encoding="utf-8")
    print("\n".join(log))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
