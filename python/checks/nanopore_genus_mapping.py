# nanopore_genus_mapping.py
# Which EMU taxa make up each Nanopore genus column of the matrix (Prusak 2026,
# code label Hanna_2025). The published run bridged the EMU names to SILVA genus
# names with a dictionary that is no longer available; this check rebuilds the
# bridge from the data alone (Preetam's point 5 of 21 Nov 2026).
#
# 120 of the 124 Nanopore genus columns have the name of an EMU genus. Each must
# equal the sum of that genus's EMU taxa in all 54 samples. For the other four
# (Escherichia-Shigella, Clostridium_sensu_stricto_12, Anaerostignum,
# Unclassified) the question is: is there exactly one grouping of the remaining
# EMU taxa that reproduces each column in all 54 samples, at the precision of
# the saved files?
#
# Method:
#   1. The 120 columns with an EMU genus name: each is checked against the sum
#      of that genus's EMU taxa, so the remaining taxa are exactly those left
#      for the four columns.
#   2. A remaining taxon can go into a column only if the column is above zero
#      in every sample where the taxon is above zero.
#   3. For each column, the candidate taxa form a 54 x n matrix A and the column
#      a vector b. A grouping is a 0/1 vector x with A x = b within tolerance.
#      - If A has full column rank and its smallest singular value s satisfies
#        s > 2 * tol * sqrt(54), no two different 0/1 vectors can both fit, so
#        a grouping that fits is the only one (a proof, not a search).
#      - Otherwise (more candidates than independent directions), any other
#        grouping must drop a set S of the members and add a set T of
#        non-members with sum(S) = sum(T) in every sample. All abundances are
#        >= 0, so T cannot be empty, and every taxon in S must be 0 wherever
#        sum(T) is 0. The search tries every non-empty T and every S that rule
#        allows.
#   4. The four columns together: a taxon goes into at most one column.
#   5. Tolerance per sample: half a unit in the last saved digit of the matrix
#      value, plus half a unit in the last saved digit of every candidate value
#      above zero in that sample, plus 1e-9 for floating-point sums. A value is
#      taken as saved to at least 15 significant digits (R's write.csv drops
#      trailing zeros, so "1" stands for 1.00000000000000).
#
# Inputs (read only):
#   - the authors' files, which are not in this repository, in one folder:
#     Nanopore_Abundance.csv and Nanopore_Taxonomy.csv (the EMU output: 286 taxa
#     in 54 samples) and "metadata_Expanding the product_Prusak.xlsx", whose
#     sheet SUCC_META maps the matrix's Sample_ID (sample_A1) to the EMU code (1A_1);
#   - the final matrix (BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv) of a corrected
#     run, by default the run of record (RUN_OF_RECORD in build_results_reference.py).
# Needs openpyxl, which is not in requirements.txt because no analysis step
# reads an Excel file (pip install openpyxl).
#
# Usage, from the repository root:
#   python python/checks/nanopore_genus_mapping.py <folder with the authors' files>
#       [--matrix-dir <corrected_<date> folder>] [--out <folder>]
# Writes into --out (default results/primer_trimmed_2026-10-06/source_checks_<today>):
#   nanopore_emu_to_matrix_genus.csv    every EMU taxon with its matrix column
#   nanopore_mapping_check_output.txt   everything printed, with the inputs' SHA-256
# Nothing is overwritten: the script stops if either file is already there.

import argparse
import csv
import datetime
import hashlib
import itertools
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import openpyxl

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "python"))
from build_results_reference import RUN_OF_RECORD, sha256_lf  # noqa: E402

MATRIX_NAME = "BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"
ABUNDANCE, TAXONOMY = "Nanopore_Abundance.csv", "Nanopore_Taxonomy.csv"
METADATA = "metadata_Expanding the product_Prusak.xlsx"
TARGETS = ["Escherichia-Shigella", "Clostridium_sensu_stricto_12", "Anaerostignum", "Unclassified"]
EMU_NAME = {"Escherichia-Shigella": "Escherichia", "Clostridium_sensu_stricto_12": "Clostridium",
            "Anaerostignum": "Anaerotignum", "Unclassified": ""}
OUT_CSV, OUT_TXT = "nanopore_emu_to_matrix_genus.csv", "nanopore_mapping_check_output.txt"

ap = argparse.ArgumentParser(description="Rebuild the EMU-to-matrix bridge of the Nanopore genus columns.")
ap.add_argument("authors_dir", type=Path, help=f"folder holding {ABUNDANCE}, {TAXONOMY} and {METADATA}")
ap.add_argument("--matrix-dir", type=Path, default=RUN_OF_RECORD, help=f"corrected_<date> folder holding {MATRIX_NAME}")
ap.add_argument("--out", type=Path, default=None,
                help="output folder (default results/primer_trimmed_2026-10-06/source_checks_<today>)")
args = ap.parse_args()
OUT = args.out or REPO / "results" / "primer_trimmed_2026-10-06" / f"source_checks_{datetime.date.today():%Y-%m-%d}"
for name in (OUT_CSV, OUT_TXT):
    if (OUT / name).exists():
        raise SystemExit(f"{OUT / name} already exists. Nothing is overwritten; give a new folder with --out.")
FINAL = args.matrix_dir / MATRIX_NAME
LOG = []


def say(line=""):
    print(line)
    LOG.append(line)


def shown(path):
    path = Path(path).resolve()
    try:
        return path.relative_to(REPO).as_posix()
    except ValueError:
        return path.name


def half_ulp(s):
    """Half a unit in the last saved digit of a number written as text, taking
    at least 15 significant digits (0 for an exact 0)."""
    s = s.strip()
    v = float(s or 0)
    if v == 0:
        return 0.0
    mant = s.lower().split("e")[0].replace("-", "").replace("+", "")
    digits = mant.replace(".", "").lstrip("0")
    nsig = max(len(digits), 15)
    e = int(np.floor(np.log10(abs(v))))
    return 0.5 * 10.0 ** (e - nsig + 1)


say("Inputs (SHA-256; text files with line endings read as LF, the workbook as saved):")
for path in (args.authors_dir / ABUNDANCE, args.authors_dir / TAXONOMY, FINAL):
    say(f"  {shown(path)}  {sha256_lf(path)}")
say(f"  {METADATA}  {hashlib.sha256((args.authors_dir / METADATA).read_bytes()).hexdigest()}")
say()

with open(args.authors_dir / TAXONOMY, encoding="utf-8-sig", newline="") as f:
    tax = list(csv.reader(f))
head = tax[0]
info = {r[0]: dict(zip(head, r)) for r in tax[1:]}
with open(args.authors_dir / ABUNDANCE, encoding="utf-8-sig", newline="") as f:
    ab = list(csv.reader(f))
codes = ab[0][1:]
raw = {r[0]: r[1:] for r in ab[1:]}
assert set(raw) == set(info), "abundance and taxonomy files list different taxa"

ws = openpyxl.load_workbook(args.authors_dir / METADATA, data_only=True, read_only=True)["SUCC_META"]
x = list(ws.iter_rows(values_only=True))
code_of = {r[x[0].index("sample")]: r[x[0].index("Sample_ID")] for r in x[1:] if r[0]}

with open(FINAL, encoding="utf-8-sig", newline="") as f:
    final = [r for r in csv.DictReader(f) if r["Paper_ID"] == "Hanna_2025"]
samples = [r["Sample_ID"] for r in final]
K = len(final)
col = [codes.index(code_of[s]) for s in samples]
assert len(set(col)) == K == 54, "sample codes do not map one to one"

gcols = [c for c in final[0] if c.startswith("g__")]
matrix_genera = [c[3:] for c in gcols if any(float(r[c] or 0) > 0 for r in final)]
genus_of = {t: info[t]["genus"] for t in info}
by_genus = defaultdict(list)
for t in info:
    by_genus[genus_of[t]].append(t)

V = {t: np.array([float(raw[t][j] or 0) for j in col]) for t in info}
U = {t: np.array([half_ulp(raw[t][j]) for j in col]) for t in info}

say(f"EMU taxa: {len(info)}; EMU genera: {len(by_genus)}; Nanopore genus columns in the matrix: {len(matrix_genera)}; samples: {K}")
say(f"EMU genera that are not matrix genus names: {sorted(g or '(no genus)' for g in by_genus if g not in matrix_genera)}")
say(f"Matrix genus columns that are not EMU genus names: {[g for g in matrix_genera if g not in by_genus]}")
zero_taxa = [t for t in info if not V[t].any()]
say(f"EMU taxa that are 0 in all 54 samples: {len(zero_taxa)}")

# 1. The columns with an EMU genus name
worst, worst_g, bad = 0.0, None, []
for g in matrix_genera:
    if g not in by_genus:
        continue
    b = np.array([float(r["g__" + g] or 0) for r in final])
    tb = np.array([half_ulp(r["g__" + g]) for r in final])
    s = sum(V[t] for t in by_genus[g])
    tol = tb + sum(U[t] for t in by_genus[g]) + 1e-9
    d = np.abs(s - b)
    if (d > tol).any():
        bad.append(g)
    if d.max() > worst:
        worst, worst_g = d.max(), g
n_verbatim = sum(1 for g in matrix_genera if g in by_genus)
say(f"\n1. Verbatim genera: {n_verbatim}; equal to the sum of their EMU taxa in all {K} samples: {n_verbatim - len(bad)}"
    f" (largest difference {worst:.3g}, {worst_g}); not equal: {bad or 'none'}")

remaining = [t for t in info if genus_of[t] not in matrix_genera]
say(f"Remaining EMU taxa (genus not a matrix genus name): {len(remaining)}, of which 0 in all samples: "
    f"{sum(1 for t in remaining if t in zero_taxa)}")
for g in sorted({genus_of[t] for t in remaining}):
    say(f"  EMU genus {g or '(no genus)'}: {sum(1 for t in remaining if genus_of[t] == g)} taxa")

target = {g: np.array([float(r["g__" + g] or 0) for r in final]) for g in TARGETS}
ttol = {g: np.array([half_ulp(r["g__" + g]) for r in final]) for g in TARGETS}

# 2. Where each remaining taxon can go. The rule is safe at this precision:
# every value above zero is far above the tolerance.
nz = min(V[t][V[t] > 0].min() for t in remaining)
say(f"Smallest value above zero among the remaining taxa: {nz:.6g}")
allowed = {t: [g for g in TARGETS if ((target[g] > 0) | (V[t] == 0)).all()] for t in remaining if t not in zero_taxa}
say("\n2. Columns each remaining taxon can go into (column > 0 wherever the taxon is > 0):")
summary = defaultdict(list)
for t, a in allowed.items():
    summary[(genus_of[t] or "(no genus)", tuple(a))].append(t)
for (g, a), ts in sorted(summary.items()):
    say(f"  EMU genus {g}: {len(ts)} taxa -> {list(a) or 'no column'}")
    if len(ts) <= 3:
        for t in ts:
            say(f"      {t}")


# 3. Column by column
def fit(members, b, tol):
    d = np.abs(sum((V[t] for t in members), np.zeros(K)) - b)
    return bool((d <= tol).all()), d


say("\n3. Each column on its own: groupings of its candidate taxa that reproduce it")
results = {}
for g in TARGETS:
    cand = [t for t in allowed if g in allowed[t]]
    b = target[g]
    tol = ttol[g] + sum((U[t] for t in cand), np.zeros(K)) + 1e-9
    A = np.column_stack([V[t] for t in cand])
    sv = np.linalg.svd(A, compute_uv=False)
    rank = int((sv > sv[0] * 1e-12).sum())
    bound = 2 * tol.max() * np.sqrt(K)
    proof = rank == len(cand) and sv[-1] > bound
    xs, *_ = np.linalg.lstsq(A, b, rcond=None)
    members = [t for t, v in zip(cand, np.round(xs).clip(0, 1)) if v == 1]
    ok, d = fit(members, b, tol)
    if not ok:  # rank-deficient: try each EMU genus among the candidates, taken whole
        for e in sorted({genus_of[t] for t in cand}):
            m = [t for t in cand if genus_of[t] == e]
            ok2, d2 = fit(m, b, tol)
            if ok2:
                members, ok, d = m, True, d2
                break
    say(f"  {g}: {len(cand)} candidate taxa, rank {rank}, smallest singular value {sv[-1]:.3g} (the proof needs > {bound:.2g})")
    say(f"      a grouping of {len(members)} taxa reproduces all {K} samples: {ok} (largest difference {d.max():.3g}; "
        f"tolerance there {tol[d.argmax()]:.3g})")
    alts, undecided = [], []
    if ok and not proof:
        # Dropping members without adding any taxon cannot fit: every member is
        # above twice the tolerance in at least one sample.
        assert all((V[t] > 2 * tol).any() for t in members), "a member is below twice the tolerance everywhere"
        non = [t for t in cand if t not in members]
        for r in range(1, len(non) + 1):
            for T in itertools.combinations(non, r):
                sT = sum(V[t] for t in T)
                zero = sT == 0
                S_cand = [t for t in members if not (V[t][zero] > 2 * tol[zero]).any()]
                if len(S_cand) > 22:
                    undecided.append(T)
                    continue
                for k in range(0, len(S_cand) + 1):
                    for S in itertools.combinations(S_cand, k):
                        alt = [t for t in members if t not in S] + list(T)
                        if fit(alt, b, tol)[0]:
                            alts.append((S, T))
        say(f"      search for other groupings: {len(alts)} found; {len(undecided)} added sets too large to search")
        for S, T in alts:
            say(f"        drop {list(S)}, add {list(T)}")
    unique_alone = ok and (proof or (not alts and not undecided))
    say(f"      the only grouping that fits, taking this column alone: {unique_alone}")
    results[g] = dict(members=members, fits=ok, unique_alone=unique_alone, alts=alts, undecided=undecided)

# 4. The four columns together: a taxon goes into at most one column
say("\n4. The four columns together")
used = defaultdict(list)
for g, r in results.items():
    for t in r["members"]:
        used[t].append(g)
twice = {t: gs for t, gs in used.items() if len(gs) > 1}
nowhere = [t for t in allowed if t not in used]
say(f"  Taxa placed in more than one column: {twice or 'none'}")
say(f"  Remaining taxa above zero that are in no column: {nowhere or 'none'}")
# An alternative for one column is impossible if it needs a taxon that another
# column's grouping, the only one for that column, already holds.
for g, r in results.items():
    held = {t for h, q in results.items() if h != g and q["unique_alone"] for t in q["members"]}
    live = [(S, T) for S, T in r["alts"] if not set(T) & held]
    live_und = [T for T in r["undecided"] if not set(T) & held]
    r["unique_joint"] = r["fits"] and not live and not live_und and not twice
    if not r["unique_alone"] and r["fits"]:
        say(f"  {g}: {len(r['alts'])} other grouping(s) taking this column alone; each needs a taxon that another column's"
            f" only grouping holds: {not live and not live_und}")
say("\nVerdict (all 54 samples, at the precision of the saved files):")
for g in TARGETS:
    r = results[g]
    whole = sorted(r["members"]) == sorted(t for t in remaining if genus_of[t] == EMU_NAME[g])
    state = ("uniquely supported" if r["unique_joint"] else "fits, but not unique") if r["fits"] else "no grouping fits"
    say(f"  {g}: {state}; = all {len(r['members'])} EMU taxa of genus {EMU_NAME[g] or '(no genus: unassigned)'}: {whole}")
    for t in r["members"]:
        say(f"      {t}")

# 5. The inferred mapping of all EMU taxa, for the supplement
OUT.mkdir(parents=True, exist_ok=True)
col_of = {t: gs[0] for t, gs in used.items()}
with open(OUT / OUT_CSV, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["emu_taxon", "emu_genus", "emu_species", "matrix_genus_column", "how_established", "samples_above_zero"])
    for t in sorted(info, key=lambda t: (genus_of[t], t)):
        if genus_of[t] in matrix_genera:
            c, how = genus_of[t], "same genus name"
        elif t in col_of:
            c = col_of[t]
            how = ("inferred from the data (only grouping that fits)" if results[c]["unique_joint"]
                   else "inferred from the data (not unique)")
        else:
            c, how = "", "in no column"
        w.writerow([t, genus_of[t], info[t].get("species", ""), c, how, int((V[t] > 0).sum())])
say(f"\nWrote {shown(OUT / OUT_CSV)} ({len(info)} EMU taxa)")
(OUT / OUT_TXT).write_text("\n".join(LOG) + "\n", encoding="utf-8")
print(f"Wrote {shown(OUT / OUT_TXT)}")
