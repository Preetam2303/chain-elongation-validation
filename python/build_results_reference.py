# build_results_reference.py
# Stage 9 of the primer-trimmed rerun: builds results_reference.json from the
# run of record in results/primer_trimmed_2026-10-06/ (the laptop's files,
# copied unchanged) and data/table_s2/. Every value is read from a results
# file: the read tracking of Stage 1, the matrices that R/04, 00 and 00b wrote,
# the tables and logs of Stage 6 (02 to 11), the figure logs (Stage 8) and
# Table S2 (12, 12b). None is typed in. Each value is rounded as the
# manuscript reports it, and each section lists the files it comes from.
#
# Each Stage 6 script's run of record is the one stage6_corrected_* folder that
# holds its log (<script>_*.log), and each figure's is the one
# figures_corrected_* folder that holds its log. The script stops if a log is
# in more than one folder.
#
# Usage, from the repository root:
#   python python/build_results_reference.py
#       writes results_reference.json from the run of record
#   python python/build_results_reference.py --check <your corrected_<date> folder>
#       reads the same values from your own run (the matrices in that folder,
#       its stage6_corrected_* and figures_corrected_* folders, and the Stage 2
#       tables in the folder above it) and lists every value that differs from
#       results_reference.json by more than half a unit in its last decimal.
#       Values whose files your run does not have are listed as not checked.

import argparse
import hashlib
import itertools
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
RUN_OF_RECORD = REPO / "results" / "primer_trimmed_2026-10-06" / "corrected_2026-10-11"
TABLE_S2_DIR = REPO / "data" / "table_s2"
REFERENCE = REPO / "results_reference.json"
sys.path.insert(0, str(REPO / "python" / "pipeline"))
sys.path.insert(0, str(REPO / "python" / "figures"))
import pipeline_settings as ps  # noqa: E402
from figure_inputs import STUDY_LABELS  # noqa: E402

TABLE_2_GENERA = ["Caproiciproducens", "Clostridium_sensu_stricto_12", "Lacticaseibacillus",
                  "Stenotrophomonas", "Lactobacillus", "Gudongella", "Acinetobacter"]


class Missing(Exception):
    """A file a value is read from is not in the run being read."""


class V(float):
    """A value read from a results file, with the number of decimals it is reported to."""

    def __new__(cls, x, decimals):
        v = super().__new__(cls, float(x))
        v.decimals = decimals
        return v


def r3(x):  # R2, SD of R2, Spearman rho, p-values, mean |SHAP|, Bray-Curtis
    return V(x, 3)


def r1(x):  # RMSE, MAE, mM C, percentages
    return V(x, 1)


def num(s):
    """A number as a log prints it, with as many decimals as it prints."""
    return V(float(s), len(s.split(".")[1]) if "." in s else 0) if re.search(r"[.eE]", s) else int(s)


def folds(cell):
    """A table's folds column: 'Study R2; Study R2; ...' (a dict) or 'R2; R2; ...' (a list)."""
    parts = [p.strip() for p in str(cell).split(";")]
    if all(" " in p for p in parts):
        return {p.rsplit(" ", 1)[0]: num(p.rsplit(" ", 1)[1]) for p in parts}
    return [num(p) for p in parts]


def one(t, /, **match):
    rows = t
    for col, value in match.items():
        rows = rows[rows[col] == value]
    if len(rows) != 1:
        raise SystemExit(f"Expected one row with {match}, found {len(rows)}.")
    return rows.iloc[0]


def found(pattern, text, what):
    """The groups of the one line of `text` that `pattern` matches."""
    hits = re.findall(pattern, text, flags=re.M)
    if len(hits) != 1:
        raise SystemExit(f"Expected one line matching {pattern!r} in {what}, found {len(hits)}.")
    return hits[0]


def section_lines(text, start, what):
    """The lines after the one line that contains `start`."""
    lines = text.split("\n")
    at = [i for i, ln in enumerate(lines) if start in ln]
    if len(at) != 1:
        raise SystemExit(f"Expected one line containing {start!r} in {what}, found {len(at)}.")
    return lines[at[0] + 1:]


def sha256_lf(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read().replace(b"\r\n", b"\n")).hexdigest()


class Run:
    """The files of one corrected run: its corrected_<date> folder, the Stage 6
    and figure folders inside it, and the run folder above it."""

    def __init__(self, corrected, table_s2=None):
        self.corrected = Path(corrected).resolve()
        self.root = self.corrected.parent
        self.table_s2 = table_s2
        self.used = []

    def file(self, *parts, base=None):
        path = (base or self.corrected).joinpath(*parts)
        if not path.exists():
            raise Missing(f"no {path}")
        self.used.append(path)
        return path

    def _folder(self, pattern, log):
        hits = [d for d in sorted(self.corrected.glob(pattern)) if list(d.glob(log))]
        if len(hits) > 1:
            raise SystemExit(f"{log} is in {len(hits)} folders: {[d.name for d in hits]}. "
                             "Keep only its run of record.")
        if not hits:
            raise Missing(f"no {pattern} folder with {log}")
        return hits[0]

    def stage6(self, script):
        return self._folder("stage6_corrected_*", f"{script}_*.log")

    def csv(self, script, name):
        return pd.read_csv(self.file(name, base=self.stage6(script)))

    def log(self, script):
        folder = self.stage6(script)
        return self.file(next(folder.glob(f"{script}_*.log")).name, base=folder).read_text(
            encoding="utf-8").replace("\r\n", "\n")

    def figure_log(self, stem):
        folder = self._folder("figures_corrected_*", f"{stem}.log")
        return self.file(f"{stem}.log", base=folder).read_text(encoding="utf-8").replace("\r\n", "\n")

    def matrix(self, name="BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv", relative=False):
        return ps.load_matrix(self.file(name), relative=relative)

    def name(self, path):
        path = Path(path).resolve()
        for base, label in ((self.root, ""), (REPO, "repository: ")):
            try:
                return label + path.relative_to(base).as_posix()
            except ValueError:
                pass
        return str(path)


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------

def primer_trimming(run):
    settings = {}
    for script in sorted((REPO / "R").glob("01?_download_dada2_*.R")):
        code = script.read_text(encoding="utf-8")
        s = {"341F": found(r'^primer_fwd <- "([ACGTNWHVRYKMSBD]+)"', code, script.name),
             "805R": found(r'^primer_rev <- "([ACGTNWHVRYKMSBD]+)"', code, script.name),
             "cutadapt_version": found(r"^# EDIT THIS: path to cutadapt \(([\d.]+),", code, script.name),
             "max_errors": int(found(r'"-e", "(\d+)", "--discard-untrimmed", "--pair-filter=any"', code, script.name)),
             "truncLen": [int(x) for x in found(r"truncLen=c\((\d+),(\d+)\)", code, script.name)]}
        settings[script.name] = s
    if len({json.dumps(s) for s in settings.values()}) != 1:
        raise SystemExit(f"The six Stage 1 scripts use different trimming settings: {settings}")
    s = next(iter(settings.values()))

    prov = pd.read_csv(REPO / "data" / "study_provenance.csv").dropna(subset=["R_script"])
    study_of = dict(zip(prov["R_script"], prov["Paper_ID"]))
    tracks = {}
    log = run.file("stage1", "dada2_overnight_log_2026-10-05.txt", base=run.root).read_text(
        encoding="utf-8").replace("\r\n", "\n")
    for script, body in re.findall(r"^===== (\S+) started .*?\n(.*?)^===== \1 FINISHED", log, flags=re.M | re.S):
        counts = re.findall(r"^\S+_1\.fastq\.gz\s+(\S+)\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+\d+\s+\d+\s*$", body, flags=re.M)
        pct = re.findall(r"^\S+_1\.fastq\.gz\s+([\d.]+)\s+[\d.]+\s+\S+\s*$", body, flags=re.M)
        if not counts or len(counts) != len(pct):
            raise SystemExit(f"Could not read the read-tracking table of {script} in the Stage 1 log.")
        tracks[study_of[script]] = pd.DataFrame({"raw": [int(c[1]) for c in counts],
                                                 "kept": [int(c[2]) for c in counts],
                                                 "pct": [float(p) for p in pct]})
    d20 = pd.read_csv(run.file("stage1", "read_tracking_Duber_2020_Master.csv", base=run.root))
    tracks[study_of["01a_download_dada2_duber2020.R"]] = d20.rename(
        columns={"primers_removed": "kept", "pct_with_primers": "pct"})
    per_study = {}
    for study in sorted(tracks):
        t = tracks[study]
        per_study[study] = {"samples": len(t), "read_pairs_in": int(t["raw"].sum()),
                            "read_pairs_with_both_primers": int(t["kept"].sum()),
                            "percent_with_both_primers": r1(100 * t["kept"].sum() / t["raw"].sum()),
                            "percent_per_sample_lowest": r1(t["pct"].min()),
                            "percent_per_sample_median": r1(t["pct"].median()),
                            "percent_per_sample_highest": r1(t["pct"].max())}

    tax = pd.read_csv(run.file("ASV_Taxonomy_Master.csv", base=run.root))
    seqs = tax["Raw_Sequence"].astype(str)
    return {
        "_note": "Both primers are cut off every read pair with cutadapt before DADA2, anchored at the start of "
                 "the read; pairs in which either primer is not found are discarded. The published run kept the "
                 "primers in the reads. The six scripts R/01a-01f use the same settings (read from them here). "
                 "The reverse primer is the one some papers call 785R.",
        "primers": {"341F": s["341F"], "805R": s["805R"]},
        "cutadapt": {"version": s["cutadapt_version"],
                     "options": f"-g ^341F -G ^805R -e {s['max_errors']} --discard-untrimmed --pair-filter=any",
                     "max_errors": s["max_errors"]},
        "dada2_truncLen_after_trimming": s["truncLen"],
        "read_pairs_by_study": per_study,
        "_read_pairs_note": "From each script's read tracking: pairs in, and pairs in which cutadapt found both "
                            "primers at the start of the reads (the only pairs DADA2 sees). Brodowski 2022 "
                            "[Ext.]'s 16 samples include the mock community, removed after DADA2.",
        "asvs_in_merged_table": len(tax),
        "asvs_starting_with_341F": int(seqs.str.contains(r"^CCTACGGG").sum()),
        "asvs_ending_with_805R_reverse_complement": int(seqs.str.contains(r"GGATTAGATACCC[CGT][AGT]GTAGTC$").sum()),
        "_asv_check_note": "The same check as R/01a-01f, on the merged and chimera-free ASV table "
                           "(ASV_Taxonomy_Master.csv): no ASV still carries a primer.",
    }


def data_construction(run):
    asv = pd.read_csv(run.file("ASV_Count_Matrix_Clean.csv", base=run.root), usecols=[0])
    asv_runs = set(asv.iloc[:, 0].astype(str).str.strip())
    n_asv = len(pd.read_csv(run.file("ASV_Count_Matrix_Clean.csv", base=run.root), nrows=0).columns) - 1
    ml = pd.read_csv(run.file("BIOTWIN_FINAL_ML_MATRIX.csv"), low_memory=False)
    ml_runs = set(ml["Sample_ID"].astype(str).str.strip())
    excluded = set(pd.read_csv(REPO / "data" / "corrections" / "excluded_samples.csv")["Sample_ID"])
    inocula = pd.read_csv(run.file("inocula_genus_profiles.csv"))
    pruned = pd.read_csv(run.file("BIOTWIN_PRUNED_ML_MATRIX.csv"), low_memory=False)
    kept = [c for c in pruned.columns if c.startswith("ASV_")]
    genus = pd.read_csv(run.file("BIOTWIN_GENUS_ML_MATRIX.csv"), nrows=0)
    s4 = pd.read_csv(run.file("stage4_genus_summary.csv"))
    s5 = pd.read_csv(run.file("stage5_genus_summary.csv"))
    gm, genus_cols = run.matrix()
    code = (REPO / "python" / "pipeline" / "00_prevalence_filter_genus.py").read_text(encoding="utf-8")
    fraction = float(found(r"^PREVALENCE = ([\d.]+)", code, "00_prevalence_filter_genus.py"))
    spike = json.loads(found(r"^SPIKE_IN_GENERA = (\[.*\])", code, "00_prevalence_filter_genus.py").replace("'", '"'))
    modelled = len(pruned)
    shared = s5[s5["in_illumina"] & s5["in_nanopore"]]["genus"]
    return {
        "stage2_asv_table": {"sequencing_runs": len(asv_runs), "asvs": n_asv,
                             "_note": "R/02 (merge, chimera removal) and R/03 (SILVA v138.1 taxonomy)."},
        "r04_matrix": {"rows": len(ml), "sequencing_runs": len(ml_runs),
                       "asv_columns": sum(c.startswith("ASV_") for c in ml.columns),
                       "runs_excluded": sorted(asv_runs & excluded),
                       "runs_without_metadata": sorted(asv_runs - ml_runs - excluded),
                       "_note": "R/04 joins the ASV table to the metadata (an inner join), applies "
                                "data/corrections/ and gives each shared inoculum run one row per reactor it "
                                "seeded. The runs without a metadata row were not in the published matrix either."},
        "inoculum_rows": {"rows": len(inocula), "sequencing_runs": inocula["Sample_ID"].nunique(),
                          "_note": "Left out of every model (data/corrections/descriptive_only_samples.csv); "
                                   "their genus profiles are in inocula_genus_profiles.csv."},
        "prevalence_filter": {"fraction": fraction, "of_modelled_illumina_samples": modelled,
                              "min_samples": math.ceil(fraction * modelled),
                              "fewest_samples_of_a_kept_asv": int((pruned[kept] > 0).sum().min()),
                              "spike_in_genera_removed_first": spike,
                              "asvs_kept": len(kept),
                              "percent_of_reads_kept": r1(s4["pct_of_all_reads"].sum()),
                              "_note": "00: the inoculum rows are left out, the spike-in genera removed, and "
                                       "an ASV is kept if it is present in at least 5% of the modelled "
                                       "Illumina samples, each counted once. Percent of reads: the kept ASVs' "
                                       "share of all reads in those samples (spike-ins removed), from "
                                       "stage4_genus_summary.csv."},
        "illumina_genera": len([c for c in genus.columns if c.startswith("g__")]),
        "nanopore_rows": int((gm["Paper_ID"] == ps.NANOPORE_ID).sum()),
        "nanopore_genera": int(s5["in_nanopore"].sum()),
        "genera_shared_across_platforms": len(shared),
        "named_genera_shared_across_platforms": int((shared != "g__Unclassified").sum()),
        "final_matrix": {"rows": len(gm), "genera": len(genus_cols), "columns": gm.shape[1],
                         "_note": "00b: the 54 Nanopore rows joined to the 112 modelled Illumina rows on genus "
                                  "columns (an outer join; a genus absent from a platform is 0)."},
    }


def table_1_sample_counts(run):
    gm, _ = run.matrix()
    inocula = pd.read_csv(run.file("inocula_genus_profiles.csv"))
    out = {}
    for study in sorted(gm["Paper_ID"].unique()):
        rows = gm[gm["Paper_ID"] == study]
        ino = inocula[inocula["Paper_ID"] == study]
        out[study] = {"modelled_rows": len(rows), "modelled_sequencing_runs": rows["Sample_ID"].nunique(),
                      "inoculum_rows": len(ino)}
    out["total"] = {"modelled_rows": len(gm), "modelled_sequencing_runs": gm["Sample_ID"].nunique(),
                    "inoculum_rows": len(inocula)}
    out["_note"] = ("A row is one sample on one reactor and day. An inoculum run that seeded several reactors "
                    "gives one inoculum row per reactor.")
    return out


def table_2_prevalence(run):
    gm, _ = run.matrix(relative=True)
    illumina = gm["Paper_ID"] != ps.NANOPORE_ID
    n_groups = gm["Paper_ID"].nunique()
    out = {"_definition": f"Detected = more than zero reads in a sample of the final {len(gm)}-row matrix; "
                          f"groups = study groups (of {n_groups}) with at least one detection. Mean relative "
                          "abundance per platform, over all of that platform's samples (Figure 4's definition)."}
    for name in TABLE_2_GENERA:
        g = gm[f"g__{name}"]
        present = g > 0
        entry = {"samples": int(present.sum()), "of": len(gm), "percent": r1(100 * present.mean()),
                 "groups": int(gm.loc[present, "Paper_ID"].nunique())}
        for platform, rows in (("illumina", illumina), ("nanopore", ~illumina)):
            entry[platform] = {"samples": int(present[rows].sum()), "of": int(rows.sum()),
                               "mean_relative_abundance_percent": V(100 * g[rows].mean(), 2)}
        out[name] = entry
    out["Clostridium_sensu_stricto_12"]["caveat"] = (
        "Not cross-platform comparable: the synonym dictionary (Methods 2.2) maps every Nanopore 'Clostridium' "
        "call onto this one clade, while the Illumina data resolve several Clostridium sensu stricto clades.")
    return out


def replicate_structure(run):
    gm, _ = run.matrix()
    units = ps.sampling_units(gm)
    key = gm["Paper_ID"] + "|" + gm["BIOREACTOR"] + "|" + gm["DAY"].astype(str)
    size = key.map(key.value_counts())
    out = {"sampling_units": int(units.nunique()), "of_rows": len(gm),
           "_note": "Same-day replicates (decision 2026-10-07): in Duber 2025 (Duber_2024) and Prusak 2026 "
                    "(Hanna_2025) a reactor's samples of one day share that day's chemistry, so they form one "
                    "sampling unit; every other row is its own unit (pipeline_settings.sampling_units)."}
    for study in ps.REPLICATE_STUDIES:
        k = key[gm["Paper_ID"] == study]
        counts = k.value_counts()
        out[study] = {"replicate_sets": len(counts), "rows": len(k),
                      "samples_per_set": sorted(set(int(c) for c in counts))}
    other = (gm["Paper_ID"].isin(ps.REPLICATE_STUDIES) == False) & (size > 1)  # noqa: E712
    out["other_rows_sharing_reactor_label_and_day"] = {
        s: int(n) for s, n in gm.loc[other, "Paper_ID"].value_counts().sort_index().items()}
    out["_other_note"] = ("Duber 2022's bottles: two bottles per trial on day 22 share a trial label but are "
                          "separate vessels with their own chemistry, so each is its own sampling unit.")
    groups = key.value_counts()
    out["rows_sharing_study_reactor_label_and_day"] = {
        "rows": int(groups[groups > 1].sum()), "groups": int((groups > 1).sum()), "of_groups": len(groups),
        "largest_group": int(groups.max()),
        "_note": "The count Methods 2.6 gave for the published matrix (124 of 184 rows; 44 of 104 reactor-days)."}
    return out


def table_3a_naive_schemes(run):
    t = run.csv("10", "10_table3a_naive_splits.csv")
    out = {}
    names = {"Single random 80/20 split (ASV level)": "single_random_80_20",
             "Random 5-fold CV (ASV level)": "random_5_fold_cv",
             "Shuffled 10-fold CV, per-fold SHAP top-25 genera": "shuffled_10_fold_cv",
             "Single-study internal 80/20 split (Duber 2024)": "single_study_internal_80_20"}
    for scheme, key in names.items():
        r = one(t, scheme=scheme)
        e = {"r2": r3(r["r2_mean"]), "rmse": r1(r["rmse_mean"]), "n": int(r["n"]),
             "features": int(r["n_features"]),
             "replicate_sets_split": int(r["replicate_sets_split"]), "replicate_sets": int(r["replicate_sets"])}
        if pd.notna(r["r2_sd"]):
            e.update(r2_sd=r3(r["r2_sd"]), rmse_sd=r1(r["rmse_sd"]), folds=folds(r["folds"]))
        if pd.notna(r["n_test"]):
            e["n_test"] = int(r["n_test"])
        out[key] = e
    out["_note"] = "10, seed 42 as in the published Table 3a. replicate_sets_split: the same-day replicate sets " \
                   "that the split puts on both sides."
    seeds = run.csv("11", "11_naive_splits_over_seeds.csv")
    over = {}
    for scheme, key in names.items():
        s = seeds[seeds["scheme"] == scheme]
        every = s[s["rows"] == "every sample"].set_index("seed")["r2"]
        per_day = s[s["rows"] == "one row per sampling day"].set_index("seed")["r2"]
        if len(every) != 20 or len(per_day) != 20:
            raise SystemExit(f"Expected 20 seeds of {scheme} in 11_naive_splits_over_seeds.csv.")
        over[key] = {"every_sample": {"median": r3(every.median()), "lowest": r3(every.min()),
                                      "highest": r3(every.max())},
                     "one_row_per_sampling_day": {"median": r3(per_day.median()), "lowest": r3(per_day.min()),
                                                  "highest": r3(per_day.max())},
                     "lower_with_one_row_per_day_in_seeds": int((per_day < every.loc[per_day.index]).sum())}
    out["seeds_0_19"] = over
    out["_seeds_note"] = "11: each naive scheme with seeds 0-19, on every sample and with one row per sampling " \
                         "day; median and range of R2."
    return out


def _loso_row(r):
    e = {"r2_mean": r3(r["r2_mean"]), "r2_sd": r3(r["r2_sd"]), "rmse": r1(r["rmse_mean"]),
         "rmse_sd": r1(r["rmse_sd"]), "n": int(r["n_rows"]), "features": int(r["n_features"])}
    e["r2_by_held_out_study"] = folds(r["folds"])
    return e


def table_3b_grouped_schemes(run):
    t2 = run.csv("02", "02_loso_tables_3b_4.csv")
    t4 = run.csv("04", "04_groupkfold_by_vessel.csv")
    t6 = run.csv("06", "06_platform_transfer.csv")
    headline = one(t2, column="full matrix", model="Tier 3: + butyrate, propionate (core set; HEADLINE)")
    out = {"loso_headline": _loso_row(headline),
           "loso_training_mean_benchmark": _loso_row(one(t2, column="full matrix", model="Training-mean benchmark")),
           "loso_shap_top25_per_fold": _loso_row(one(t2, model="LOSO, core set, per-fold SHAP top-25 genera"))}
    out["loso_headline"]["_note"] = ("The headline: leave-one-study-out on the full matrix with the core set "
                                     "(Table 4's full-matrix tier 3 is the same number).")
    for scheme, key in (("Mode 1: CSTR Only (B1, B2 Vessels)", "vessel_grouped_cstr_only"),
                        ("Mode 2: All Vessels (CSTR + Batch R1-R9)", "vessel_grouped_all_vessels")):
        for model, suffix in (("XGBoost, core set + per-fold SHAP top-25 genera", ""),
                              ("Training-mean benchmark", "_training_mean_benchmark")):
            r = one(t4, scheme=scheme, model=model)
            out[key + suffix] = {"r2_mean": r3(r["r2_mean"]), "r2_sd": r3(r["r2_sd"]), "rmse": r1(r["rmse_mean"]),
                                 "rmse_sd": r1(r["rmse_sd"]), "n": int(r["n_rows"]), "vessels": int(r["n_vessels"])}
    out["_vessel_note"] = ("04: GroupKFold by physical reactor (5 folds); ties between vessels with the same number "
                           "of rows are broken by name, so the folds are the same on every machine.")
    for scheme, key in (("Platform transfer, raw counts", "platform_transfer_raw_counts"),
                        ("Platform transfer, relative abundance", "platform_transfer_relative_abundance"),
                        ("Platform transfer", "platform_transfer_training_mean_benchmark")):
        r = one(t6, scheme=scheme)
        e = {"r2": r3(r["r2"]), "rmse": r1(r["rmse"]), "train_n": int(r["train_n"]), "test_n": int(r["test_n"])}
        if pd.notna(r["spearman_rho_day_means"]):
            e["spearman_rho_sampling_day_means"] = r3(r["spearman_rho_day_means"])
        out[key] = e
    return out


def table_4_substrate_ladder(run):
    t = run.csv("02", "02_loso_tables_3b_4.csv")
    out = {"_note": "02: leave-one-study-out, adding predictor tiers, on the full matrix (every study) and on "
                    "the Illumina studies alone. Only the Illumina-only column adds the compounds measured in the "
                    "Illumina studies alone (its tiers 3 and 4)."}
    for column, key in (("full matrix", "full_matrix"), ("Illumina-only", "illumina_only")):
        rows = t[(t["column"] == column) & t["model"].str.startswith(("Tier", "Training-mean"))]
        out[key] = {}
        for _, r in rows.iterrows():
            name = "training_mean_benchmark" if r["model"] == "Training-mean benchmark" else \
                "tier_" + r["model"].split(":")[0].split()[1]
            out[key][name] = {"description": r["model"], **_loso_row(r)}
    out["chemistry_only_no_genera"] = _loso_row(one(t, model="Chemistry only: core set without genera"))
    return out


def table_5_intrastudy_reactor_transfer(run):
    t = run.csv("05", "05_table5_reactor_pairs.csv")
    gm, _ = run.matrix()
    out = {"_note": "Train on one reactor, test on the other reactor of the same study (pairs named in "
                    "pipeline_settings.REACTOR_PAIRS). The two reactors of every pair ran under different "
                    "conditions by design (conditions below), so these are not matched-condition transfers. "
                    "Spearman rho on sampling-unit means."}
    for _, r in t.iterrows():
        study = r["Paper"]
        e = {"train_reactor": r["Train_Reactor"], "test_reactor": r["Test_Reactor"], "train_n": int(r["Train_n"]),
             "test_n": int(r["Test_n"]), "test_sampling_units": int(r["Test_sampling_units"]),
             "r2": r3(r["R2"]), "rmse": r1(r["RMSE (mM C)"]), "spearman_rho": r3(r["Spearman_rho_unit_means"])}
        for reactor in (r["Train_Reactor"], r["Test_Reactor"]):
            rows = gm[(gm["Paper_ID"] == study) & (gm["BIOREACTOR"] == reactor)]
            cond = {"samples": len(rows),
                    "feed": sorted(set(rows["Feed_Complexity"].astype(str) + " / " +
                                       rows["Primary_Carbon_Signature"].astype(str)))}
            for col in ("PH", "TEMP", "HRT"):
                cond[col] = {f"{v:g}": int(n) for v, n in rows[col].value_counts().sort_index().items()}
            e[f"conditions_{reactor}"] = cond
        out[study] = e
    lines = section_lines(run.log("05"), "Not in Table 5:", "05's log")
    out["not_in_table_5"] = dict(re.findall(r"^  ([A-Z][A-Za-z]+_\d{4}\w*): (.+)$", "\n".join(lines), flags=re.M))
    return out


def figure_3_leakage_tightening_progression(run):
    d = one(run.csv("07d", "07d_figure3_stage1.csv"), Paper="Duber_2024")
    c = one(run.csv("07c", "07c_figure3_stage2.csv"), Paper="Duber_2024")
    s3 = one(run.csv("05", "05_table5_reactor_pairs.csv"), Paper="Duber_2024")
    top = section_lines(run.log("07c"), "TOP 10 LOCAL FEATURES", "07c's log")
    shap = re.findall(r"^\s*(\S.*?)\s+(\d+\.\d+)\s*$", "\n".join(top[1:11]), flags=re.M)
    if len(shap) != 10:
        raise SystemExit("Could not read 07c's top 10 SHAP features.")
    return {
        "_note": "The same Duber 2025 (Duber_2024) split, B1 (30 samples) -> B2 (21), evaluated three times as "
                 "leakage controls are tightened. No sample changes between stages.",
        "stage_1_unscaled_random_forest": {
            "r2": V(d["R2"], 8), "rmsd": V(d["RMSD"], 3), "train_n": int(d["Train_Size"]),
            "test_n": int(d["Test_Size"]),
            "_note": "07d, scikit-learn random forest (seeded). Its R2 rounds to 0.155 at three decimals but sits "
                     "next to the rounding boundary, so it is given to eight."},
        "stage_2_intermediate_unscaled_xgboost": {
            "r2": r3(c["R2"]), "rmse": r1(c["RMSE"]), "train_n": int(c["Train_Size"]), "test_n": int(c["Test_Size"]),
            "top_10_mean_abs_shap": {f: num(v) for f, v in shap},
            "_note": "07c: no scaling at all, and its broader predictor list keeps DAY and the alcohols."},
        "stage_3_fully_leakage_controlled": {
            "r2": r3(s3["R2"]), "rmse": r1(s3["RMSE (mM C)"]), "train_n": int(s3["Train_n"]),
            "test_n": int(s3["Test_n"]),
            "_note": "05: the core set, predictors scaled on the training rows only (Table 5's Duber 2025 row)."},
    }


def permutation_control(run):
    t = run.csv("03", "03_permutation_tests.csv")
    by = run.csv("03", "03_permutation_tests_by_study.csv")
    out = {"_note": "03: the headline LOSO refitted 999 times for each shuffle of caproate by sampling day. "
                    "k = shuffles that do at least as well as the true labels (R2 as high, or RMSE as low); "
                    "p = (k + 1) / 1000, so 0.001 is the smallest possible. Shuffle i is seeded with i, so both "
                    "null distributions reproduce exactly; the Monte Carlo standard error is sqrt(p (1 - p) / 999)."}
    for shuffle, key in (("across studies", "across_studies"), ("within each study", "within_each_study")):
        r = one(t, shuffle=shuffle)
        n = int(r["n_shuffles"])
        e = {"shuffles": n, "observed_r2": r3(r["observed_r2"]), "observed_rmse": r1(r["observed_rmse"]),
             "training_mean_benchmark_rmse": r1(r["training_mean_benchmark_rmse"]),
             "null_r2_mean": r3(r["null_r2_mean"]), "null_r2_sd": r3(r["null_r2_sd"]),
             "null_r2_median": r3(r["null_r2_median"]), "k_r2": int(r["k_r2"]), "p_r2": r3(r["p_r2"]),
             "p_r2_monte_carlo_se": r3(math.sqrt(r["p_r2"] * (1 - r["p_r2"]) / n)),
             "null_rmse_mean": r1(r["null_rmse_mean"]), "null_rmse_sd": r1(r["null_rmse_sd"]),
             "k_rmse": int(r["k_rmse"]), "p_rmse": r3(r["p_rmse"]),
             "p_rmse_monte_carlo_se": r3(math.sqrt(r["p_rmse"] * (1 - r["p_rmse"]) / n))}
        e["by_held_out_study"] = {
            s["held_out_study"]: {"observed_r2": r3(s["observed_r2"]), "observed_rmse": r1(s["observed_rmse"]),
                                  "k_r2": int(s["k_r2"]), "p_r2": r3(s["p_r2"]),
                                  "k_rmse": int(s["k_rmse"]), "p_rmse": r3(s["p_rmse"])}
            for _, s in by[by["shuffle"] == shuffle].iterrows()}
        out[key] = e
    out["across_studies"]["_note"] = "Caproate shuffled by sampling day across all studies, as in the published test."
    out["within_each_study"]["_note"] = ("Caproate shuffled by sampling day within each study, so each study keeps "
                                         "its own values (added 2026-10-11).")
    return out


def algorithm_selection(run):
    t = run.csv("07b", "07b_algorithm_comparison.csv")
    seeds = run.csv("07b", "07b_algorithm_comparison_over_seeds.csv")
    out = {"_note": "07b: XGBoost and LightGBM with the same predictors on three schemes (seed 42), and the two "
                    "naive schemes with seeds 0-19."}
    for scheme in t["scheme"].unique():
        e = {}
        for _, r in t[t["scheme"] == scheme].iterrows():
            m = {"r2": r3(r["r2_mean"]), "rmse": r1(r["rmse_mean"]), "n": int(r["n_rows"]),
                 "features": int(r["n_features"])}
            if pd.notna(r["r2_sd"]):
                m.update(r2_sd=r3(r["r2_sd"]), rmse_sd=r1(r["rmse_sd"]))
            if str(r["folds"]) != "nan" and r["scheme"].startswith("LOSO"):
                m["r2_by_held_out_study"] = folds(r["folds"])
            e[r["model"]] = m
        out[scheme] = e
    over = {}
    for scheme in seeds["scheme"].unique():
        s = seeds[seeds["scheme"] == scheme].pivot(index="seed", columns="model", values="r2")
        if len(s) != 20:
            raise SystemExit(f"Expected 20 seeds of {scheme} in 07b_algorithm_comparison_over_seeds.csv.")
        over[scheme] = {m: {"median": r3(s[m].median()), "lowest": r3(s[m].min()), "highest": r3(s[m].max())}
                        for m in ("XGBoost", "LightGBM")}
        over[scheme]["lightgbm_higher_in_seeds"] = int((s["LightGBM"] > s["XGBoost"]).sum())
        over[scheme]["median_difference_lightgbm_minus_xgboost"] = r3((s["LightGBM"] - s["XGBoost"]).median())
    out["seeds_0_19"] = over
    log = run.log("07a")
    scopes = re.findall(r"^(\d)\. (.+?):\n\s+Features: (\d+)\n\s+R-squared \(R2\): (-?\d+\.\d+)\n\s+RMSE: "
                        r"(\d+\.\d+) mM C", log, flags=re.M)
    if len(scopes) != 3:
        raise SystemExit("Could not read 07a's three models.")
    top = section_lines(log, "Top 10 features by mean |SHAP|", "07a's log")
    out["lightgbm_random_forest_07a"] = {
        "_note": "07a, LightGBM in random-forest mode (bagging, seeded) on the single 80/20 split at ASV level, "
                 "with three predictor scopes. The Maximum scope adds the five compounds measured only in the "
                 "Illumina studies; its attribution is led by caprylate (the target-leakage demonstration of "
                 "Methods 2.7).",
        **{name: {"features": int(f), "r2": num(r2), "rmse": num(rmse)} for _, name, f, r2, rmse in scopes},
        "top_10_features_of_the_best_model": re.findall(r"^\s+\d+\. (\S+)$", "\n".join(top[:10]), flags=re.M)}
    return out


def ecological_robustness_claim(run):
    log = run.log("09")
    n_folds, n_samples, n_genera = found(r"--\s+(\d+) folds on the (\d+)-sample, (\d+)-genus matrix", log, "09's log")
    loso, gkf, s10 = found(r"LOSO = (\d+) folds \| GroupKFold-by-vessel = (\d+) \| Shuffled 10-fold = (\d+)", log,
                           "09's log")
    ranks = run.csv("09", "guild_membership_shap_ranks.csv").set_index("genus")
    n = int(n_folds)

    def genus(g):
        r = ranks.loc[g]
        return {"median_rank": num(f"{r['median_rank']:g}"), "best_rank": int(r["best"]),
                "worst_rank": int(r["worst"]), "top_5_folds": int(r["top5"]), "top_10_folds": int(r["top10"]),
                "top_25_folds": int(r["top25"]), "top_10_in_loso_folds": int(r["loso_top10"]),
                "top_10_in_vessel_folds": int(r["gkf_top10"]), "top_10_in_shuffled_10_fold": int(r["s10_top10"]),
                "of_folds": n, "zero_attribution_without_nanopore": bool(r["zero_when_no_nano"]),
                "shap_direction": r3(r["shap_dir"]), "raw_spearman_rho": r3(r["raw_rho"])}

    rank2 = re.findall(r"^\s+(g__.+?)\s+rank 2 in (\d+)/(\d+) folds$", log, flags=re.M)
    cross = section_lines(log, "CROSS-CHECK: Illumina-only", "09's log")
    within = run.csv("09", "guild_producer_within_study_correlation.csv")
    return {
        "_definition": f"09: each genus's rank by mean |SHAP| in each of {n} folds ({loso} LOSO, {gkf} "
                       f"GroupKFold-by-vessel, {s10} shuffled 10-fold) on the {n_samples}-sample, "
                       f"{n_genera}-genus matrix.",
        "folds": {"total": n, "loso": int(loso), "vessel": int(gkf), "shuffled_10_fold": int(s10)},
        "Caproiciproducens": genus("g__Caproiciproducens"),
        "Clostridium_sensu_stricto_12": genus("g__Clostridium_sensu_stricto_12"),
        "top_10_in_every_fold": sorted(ranks.index[ranks["top10"] == n]),
        "top_10_in_at_least_80_percent_of_folds": sorted(ranks.index[ranks["top10"] >= 0.8 * n]),
        "rank_2_folds": {g: int(k) for g, k, _ in rank2},
        "illumina_only_109_genus_model_top_10": re.findall(r"^\s+\d+\. (g__.+)$", "\n".join(cross[1:12]), flags=re.M),
        "caproiciproducens_vs_caproate_within_study": {
            r["study"]: {"rows": int(r["n_rows"]), "sampling_units": int(r["n_sampling_units"]),
                         "mean_relative_abundance_percent": r1(r["mean_rel_abundance_pct"]),
                         "spearman_rho": r3(r["spearman_rho"]), "p": V(r["p_value"], 3 if r["p_value"] > 0.001 else 14)}
            for _, r in within.iterrows()},
        "_within_note": "Spearman rho of Caproiciproducens relative abundance with measured caproate inside each "
                        "study, on sampling-unit means.",
    }


def shared_subspace_transfer_control(run):
    log = run.log("08")
    space = {}
    for label, key in (("Total genus columns in matrix", "genera_total"), ("Detected on Illumina", "detected_illumina"),
                       ("Detected on Nanopore", "detected_nanopore"), ("SHARED \\(detected on both\\)", "shared"),
                       ("Illumina-only \\(forced to 0 at test\\)", "illumina_only"),
                       ("Nanopore-only \\(no training signal\\)", "nanopore_only")):
        space[key] = int(found(rf"^\s+{label}\s+: (\d+)$", log, "08's log"))
    space["shared_mass_illumina_percent"] = num(found(r"^\s+Illumina samples :\s+([\d.]+)%$", log, "08's log"))
    space["shared_mass_nanopore_percent"] = num(found(r"^\s+Nanopore samples :\s+([\d.]+)%$", log, "08's log"))
    unc = found(r"^\s+Illumina :\s+([\d.]+)%\s+Nanopore :\s+([\d.]+)%$", log, "08's log")
    space["unclassified_mass_percent"] = {"illumina": num(unc[0]), "nanopore": num(unc[1])}
    listed = section_lines(log, "shared genera:", "08's log")
    space["shared_genera"] = [ln.strip() for ln in itertools.takewhile(lambda ln: ln.startswith("    g__"), listed)]
    if len(space["shared_genera"]) != space["shared"]:
        raise SystemExit("08's log lists a different number of shared genera than it counts.")

    def table(start, pattern, fields):
        rows = {}
        lines = section_lines(log, start, "08's log")
        for ln in lines[1:] if lines and lines[0].startswith("=====") else lines:
            if ln.startswith("====="):
                break
            m = re.match(pattern, ln)
            if m:
                label = re.sub(r"\s*(<-- KEY|\[BASELINE\])", "", m.group(1))
                label = re.sub(r"\s+", " ", re.sub(r"^(S\d)\s+", r"\1: ", label)).strip()
                rows[label] = {f: num(x) for f, x in zip(fields, m.groups()[1:])}
        if not rows:
            raise SystemExit(f"Could not read the table after {start!r} in 08's log.")
        return rows

    transfer = r"^\s+((?:S\d\s+)?.+?)\s+(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
    guard = r"^\s+(.+?)\s+(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+\["
    cols = ("features", "r2", "rmse", "mae", "spearman_rho")
    return {
        "_note": "08 (Table S6): Illumina -> Nanopore transfer restricted to feature scopes, with two guards. "
                 "Genera as relative abundance over all 205 genera, then subset (or renormalized within the subset).",
        "shared_feature_space": space,
        "illumina_to_nanopore": table("[global rel. abundance over all 205, then subset]", transfer, cols),
        "illumina_to_nanopore_renormalized": table("[rel. abundance RENORMALIZED within retained subset]",
                                                   transfer, cols),
        "guard_within_illumina_loso": table("GUARD (a): WITHIN-ILLUMINA LOSO", guard, ("features", "r2_mean", "r2_sd")),
        "guard_reverse_nanopore_to_illumina": table("GUARD (b): REVERSE TRANSFER", transfer, cols),
        "consistency_with_06": found(r"^\s+06 gave \S+\s+->\s+(\w+)$", log, "08's log"),
    }


def replicate_checks(run):
    bc = run.csv("11", "11_replicate_bray_curtis.csv")
    per_day = run.csv("11", "11_one_row_per_day.csv")
    out = {"_note": "11: how similar same-day replicates are, and the headline and Table 3a with one row per "
                    "sampling day (each day's replicates averaged) instead of every sample.",
           "bray_curtis": {}}
    for study in bc["study"].unique():
        rows = bc[bc["study"] == study]
        out["bray_curtis"][study] = {
            r["pairs"]: {"pairs": int(r["n_pairs"]), "median": r3(r["median"]), "q25": r3(r["q25"]),
                         "q75": r3(r["q75"]), "lowest": r3(r["min"]), "highest": r3(r["max"])}
            for _, r in rows.iterrows()}
        out["bray_curtis"][study]["samples_whose_most_similar_sample_is_a_same_day_replicate"] = int(
            rows["samples_nearest_to_own_replicate"].iloc[0])
        out["bray_curtis"][study]["samples"] = int(rows["samples"].iloc[0])
    out["one_row_per_sampling_day"] = {}
    for (scheme, model), rows in per_day.groupby(["scheme", "model"], sort=False):
        e = {}
        for _, r in rows.iterrows():
            m = {"r2": r3(r["r2_mean"]), "rmse": r1(r["rmse_mean"]), "n": int(r["n"])}
            if pd.notna(r["r2_sd"]):
                m["r2_sd"] = r3(r["r2_sd"])
            e[r["rows"]] = m
        out["one_row_per_sampling_day"][scheme if model == "XGBoost" else f"{scheme}: {model}"] = e
    return out


def figures(run):
    out = {"_note": "Values as the figure logs print them (Stage 8, run_figures.py); every input of each figure, "
                    "with its SHA-256, is listed in its log."}
    f4 = run.figure_log("figure4_cohort_comparison")
    cohorts = found(r"^Cohort\s+Illumina \(n=(\d+)\)\s+Nanopore \(n=(\d+)\)\s*$", f4, "Figure 4's log")
    x = r"\s+(-?\d+\.\d+)"
    means = re.findall(rf"^(\S.*?){x * 6}\s*$", f4, flags=re.M)
    out["figure_4_cohorts"] = {
        "samples": {"illumina": int(cohorts[0]), "nanopore": int(cohorts[1])},
        **{v: {"illumina": {"n": int(float(a)), "mean": num(b), "median": num(c)},
               "nanopore": {"n": int(float(d)), "mean": num(e), "median": num(f)}}
           for v, a, b, c, d, e, f in means}}
    f5 = run.figure_log("figure5_shap_genus")
    top15 = re.findall(r"^\s*(g__.+?)\s+(\d+\.\d+)\s+(\d+)\.000\s+(\d+)\.000\s*$", f5, flags=re.M)
    shown = found(r"^\[Figure 5\] Shown: .*: (g__.+)$", f5, "Figure 5's log")
    out["figure_5_top_10_genera"] = {
        "_note": "One XGBoost fit on all 166 modelled samples, TreeSHAP on the same samples; the top 10 genera by "
                 "mean |SHAP| (mM C), the Unclassified bin left out.",
        "shown": [g.strip() for g in shown.split(",")],
        "top_15_genera_mean_abs_shap": {g: {"mean_abs_shap": num(s), "rank_among_all_predictors": int(a),
                                            "rank_among_genera": int(b)} for g, s, a, b in top15}}
    f6 = run.figure_log("figure6_attribution_vs_spread")
    plotted, n_folds, n_groups = found(r"^\[Figure 6\] (\d+) genera plotted \((\d+) folds, (\d+) study groups\)", f6,
                                       "Figure 6's log")
    named = re.findall(r"^\s*(g__.+?)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(True|False)\s*$", f6, flags=re.M)
    hi = found(r"median rank <=10: mean groups present = ([\d.]+) \(n=(\d+)\)", f6, "Figure 6's log")
    lo = found(r"median rank  >25: mean groups present = ([\d.]+) \(n=(\d+)\)", f6, "Figure 6's log")
    out["figure_6_attribution_vs_spread"] = {
        "genera_plotted": int(plotted), "folds": int(n_folds), "study_groups": int(n_groups),
        "named_genera": {g: {"groups": int(a), "samples": int(b), "median_rank": num(m), "shap_direction": num(d),
                             "zero_attribution_without_nanopore": z == "True"}
                         for g, a, b, m, d, z in sorted(named, key=lambda x: (float(x[3]), x[0]))},
        "mean_groups_of_genera_with_median_rank_at_most_10": {"mean": num(hi[0]), "genera": int(hi[1])},
        "mean_groups_of_genera_with_median_rank_above_25": {"mean": num(lo[0]), "genera": int(lo[1])}}
    s5 = run.figure_log("figure_s5_combined_shap")
    top = re.findall(r"^\s*(\S.*?)\s+(\d+\.\d+)\s+(\d+)\.000\s*$", s5, flags=re.M)
    out["figure_s5_top_15_predictors"] = {f: {"mean_abs_shap": num(s), "rank": int(r)} for f, s, r in top}
    for key, n in (("figure_5_top_10_genera", 10), ("figure_s5_top_15_predictors", 15)):
        got = out[key]["shown"] if key.startswith("figure_5") else out[key]
        if len(got) != n:
            raise SystemExit(f"Read {len(got)} entries for {key}, expected {n}.")
    return out


def table_s2(run):
    if run.table_s2 is None:
        raise Missing("no Table S2 folder")
    ranking = pd.read_csv(run.file("12_asv_shap_ranking.csv", base=run.table_s2))
    tables = sorted(run.table_s2.glob("TableS2_*.csv"))
    if len(tables) != 1:
        raise Missing(f"{len(tables)} TableS2_*.csv files in {run.table_s2}")
    t = pd.read_csv(run.file(tables[0].name, base=run.table_s2))
    capro = t[t["included_as"].str.contains("Caproiciproducens 90%")].sort_values("cum_pct_of_genus_reads")
    lacta = capro["closest_relatives"].str.contains("lactatifermentans") & capro["at_least_99"]
    other = capro["at_least_99"] & ~lacta
    below = ~capro["at_least_99"]
    top20 = t[t["included_as"].str.startswith("SHAP top 20")].sort_values("rank")

    def share(rows):
        return r1(capro.loc[rows, "pct_of_genus_reads"].sum())

    return {
        "_note": "12 ranks every kept ASV by mean |SHAP| in one XGBoost fit on all 112 modelled Illumina samples; "
                 "12b builds the table from the BLAST of its FASTA (RID below). A top hit of 99% or more is the "
                 "closest known relative, not a species identification.",
        "asvs_ranked": len(ranking),
        "blast_rid": str(t["RID"].iloc[0]),
        "top_20": {r["ASV_ID"]: {"rank": int(r["rank"]), "silva_taxon": r["silva_taxon"],
                                 "mean_abs_shap": r3(r["mean_abs_shap"]), "studies": int(r["n_studies"]),
                                 "closest_relatives": r["closest_relatives"], "identity_pct": V(r["identity_pct"], 2),
                                 "query_cover_pct": int(r["query_cover_pct"])}
                   for _, r in top20.iterrows()},
        "caproiciproducens_asvs_with_90_percent_of_reads": {
            "asvs": len(capro), "percent_of_genus_reads": r1(capro["cum_pct_of_genus_reads"].max()),
            "at_least_99_to_C_lactatifermentans": {"asvs": sorted(capro.loc[lacta, "ASV_ID"]),
                                                   "percent_of_genus_reads": share(lacta)},
            "at_least_99_to_another_species": {"asvs": int(other.sum()), "percent_of_genus_reads": share(other)},
            "below_99_to_any_type_strain": {"asvs": int(below.sum()), "percent_of_genus_reads": share(below)}},
    }


SECTIONS = [primer_trimming, data_construction, table_1_sample_counts, table_2_prevalence, replicate_structure,
            table_3a_naive_schemes, table_3b_grouped_schemes, table_4_substrate_ladder,
            table_5_intrastudy_reactor_transfer, figure_3_leakage_tightening_progression, permutation_control,
            algorithm_selection, ecological_robustness_claim, shared_subspace_transfer_control, replicate_checks,
            figures, table_s2]


def about(run):
    versions = (run.stage6("02") / "versions.txt").read_text(encoding="utf-8").replace("\r\n", "\n")
    python = found(r"^Python (\S+) on (\S+)$", versions, "versions.txt")
    stage6 = {}
    for script in ("02", "03", "04", "05", "06", "07a", "07b", "07c", "07d", "08", "09", "10", "11"):
        stage6[script] = run.name(run.stage6(script))
    matrices = {p.name: sha256_lf(p) for p in sorted(run.corrected.glob("*.csv"))}
    return {
        "purpose": "The corrected, primer-trimmed run's values for the manuscript, with the results files each "
                   "comes from and how exactly it reproduces. Built by python/build_results_reference.py from "
                   "the run of record; no value is typed in by hand.",
        "manuscript": "Quantitative Fragility and Ecological Robustness in Machine-Learning-Predicted Chain "
                      "Elongation: A Multi-Study, Multi-Platform Validation Framework",
        "run_of_record": {
            "folder": "results/primer_trimmed_2026-10-06 (copied unchanged from the laptop run; see its README)",
            "matrices": f"{run.name(run.corrected)} (R/04_metadata_merge.R, 00 and 00b, 2026-10-11)",
            "matrix_sha256_lf": matrices,
            "stage6_folder_of_each_script": stage6,
            "figures": run.name(run._folder("figures_corrected_*", "figure1_*.log")),
            "table_s2": "data/table_s2 (12 and 12b, 2026-10-12)",
        },
        "reference_environment": {"python": python[0], "os": python[1],
                                  "packages": "exactly as pinned in requirements.txt"},
        "study_names": {k: v for k, v in STUDY_LABELS.items()},
        "precision": "Each value is rounded as the manuscript reports it: R2, Spearman rho, p-values, mean |SHAP| "
                     "and Bray-Curtis to 3 decimals; RMSE, MAE, concentrations (mM C) and percentages to 1; counts "
                     "exactly. Values read from a log keep the decimals the log prints.",
        "tolerance": {
            "exact": "Every model of the corrected pipeline is deterministic. XGBoost runs without row or column "
                     "subsampling; the few random parts are seeded (07a's LightGBM random-forest mode, 07d's "
                     "random forest, the 999 shuffles of 03, the seeds 0-19 of 07b and 11), and the vessel folds "
                     "of 04 and 09 break ties by name. With the pinned versions every value reproduces to every "
                     "decimal given, on any machine: --check on the cloud's own runs of the same scripts (Linux, "
                     "from its own copy of the matrices) finds no value that differs from the laptop run of "
                     "record (Windows 11).",
            "monte_carlo": "The p-values of 03 reproduce exactly, because every shuffle is seeded, but they are "
                           "estimates over 999 shuffles; their Monte Carlo standard error is given with them.",
            "check": "python python/build_results_reference.py --check <your corrected_<date> folder> lists every "
                     "value of your run that differs from this file by more than half a unit in its last decimal.",
        },
        "earlier_values": "The published run's values (before primer trimming, with the earlier settings and "
                          "predictors) are in this file as it stood at tag draft-v1-2026-10-04 (commit 29635f7). "
                          "They reproduce with the code at commit 844f478.",
        "sources_note": "Paths in each section's _source are relative to results/primer_trimmed_2026-10-06 unless "
                        "marked as repository paths.",
    }


def build(run, strict):
    out = {}
    for fn in SECTIONS:
        run.used = []
        try:
            values = fn(run)
        except Missing as e:
            if strict:
                raise SystemExit(f"{fn.__name__}: {e}")
            out[fn.__name__] = {"_not_checked": str(e)}
            continue
        out[fn.__name__] = {"_source": sorted({run.name(p) for p in run.used}), **values}
    return {"_about": about(run) if strict else {}, **out}


def plain(x):
    """The tree as JSON values: each V rounded to its decimals."""
    if isinstance(x, V):
        return round(float(x), x.decimals) if x.decimals else int(round(float(x)))
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, dict):  # notes and sources first
        return {k: plain(x[k]) for k in sorted(x, key=lambda k: not str(k).startswith("_"))}
    if isinstance(x, (list, tuple)):
        return [plain(v) for v in x]
    return x


def compare(ref, new, path, report):
    if isinstance(ref, dict):
        if not isinstance(new, dict):
            report["differ"].append(f"{path}: {new!r} instead of a table")
            return
        if "_not_checked" in new:
            report["not_checked"].append(f"{path}: {new['_not_checked']}")
            return
        for k, v in ref.items():
            if k.startswith("_"):
                continue
            if k not in new:
                report["differ"].append(f"{path}/{k}: missing from your run")
            else:
                compare(v, new[k], f"{path}/{k}", report)
        for k in new:
            if not k.startswith("_") and k not in ref:
                report["differ"].append(f"{path}/{k}: in your run, not in results_reference.json")
        return
    if isinstance(ref, list):
        if not isinstance(new, (list, tuple)) or len(new) != len(ref):
            report["differ"].append(f"{path}: {plain(new)!r} instead of {ref!r}")
            return
        for i, (a, b) in enumerate(zip(ref, new)):
            compare(a, b, f"{path}[{i}]", report)
        return
    report["values"] += 1
    if isinstance(new, V) and not isinstance(ref, (bool, str)):
        if abs(float(new) - ref) > 0.5 * 10 ** -new.decimals + 1e-9:
            report["differ"].append(f"{path}: {round(float(new), new.decimals + 3)} in your run, {ref} here")
    elif plain(new) != ref:
        report["differ"].append(f"{path}: {plain(new)!r} in your run, {ref!r} here")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", metavar="CORRECTED_DIR",
                        help="compare your own run's corrected_<date> folder with results_reference.json")
    args = parser.parse_args()
    if not args.check:
        run = Run(RUN_OF_RECORD, table_s2=TABLE_S2_DIR)
        tree = plain(build(run, strict=True))
        with open(REFERENCE, "w", encoding="utf-8", newline="\n") as f:
            json.dump(tree, f, indent=1, ensure_ascii=False)
            f.write("\n")
        n = sum(1 for _ in re.finditer(r'": (?:-?\d|true|false|")', json.dumps(tree)))
        print(f"Wrote {REFERENCE} ({len(tree) - 1} sections, about {n} values).")
        return
    corrected = Path(args.check)
    if not corrected.is_dir():
        raise SystemExit(f"{corrected} is not a folder.")
    s2 = [d for d in sorted(corrected.glob("stage6_corrected_*")) if list(d.glob("12_*.log"))]
    run = Run(corrected, table_s2=s2[-1] if s2 else None)
    with open(REFERENCE, encoding="utf-8") as f:
        ref = json.load(f)
    new = build(run, strict=False)
    total = {"values": 0, "differ": [], "not_checked": []}
    for name, values in ref.items():
        if name.startswith("_"):
            continue
        report = {"values": 0, "differ": [], "not_checked": []}
        compare(values, new.get(name, {"_not_checked": "section not built"}), name, report)
        state = "not checked" if report["not_checked"] and not report["values"] else \
            ("OK" if not report["differ"] else f"{len(report['differ'])} differ")
        print(f"{name}: {report['values']} values compared, {state}")
        for line in report["differ"] + [f"not checked: {x}" for x in report["not_checked"]]:
            print(f"  {line}")
        for k in total:
            total[k] += report[k]
    print(f"\n{total['values']} values compared; {len(total['differ'])} differ; "
          f"{len(total['not_checked'])} sections or tables not checked.")
    sys.exit(1 if total["differ"] else 0)


if __name__ == "__main__":
    main()
