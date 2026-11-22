# trace_caproate.py
# One caproate value per study, traced from the authors' own file to the final
# matrix of a corrected run, so the units can be reported at the level verified
# (Preetam's point 6 of 21 Nov 2026). The matrix is in mM C. Five of the seven
# source columns traced here are already in mM C; the Duber 2022 and Prusak 2026
# columns are in g/L and are converted with
#     mM C = g/L x 1000 / (M / nC),   M = 116.16 g/mol and nC = 6 for caproate,
# the formula used for every study. A value agrees if it is within 0.51 mM C of
# the source (the three Brodowski studies are whole numbers in the matrix).
# Only caproate is traced: this says nothing about the units of other metabolite
# columns. The source rows are the ones the audit of 6 Oct 2026 mapped from each
# matrix row to its row in the author's file.
#
# Inputs (read only):
#   - the authors' files, which are not in this repository, in one folder: the
#     seven workbooks named in TRACES below;
#   - the final matrix (BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv) of a corrected
#     run, by default the run of record (RUN_OF_RECORD in build_results_reference.py).
# Needs openpyxl, which is not in requirements.txt because no analysis step
# reads an Excel file (pip install openpyxl).
#
# Usage, from the repository root:
#   python python/checks/trace_caproate.py <folder with the authors' files>
#       [--matrix-dir <corrected_<date> folder>] [--out <folder>]
# Writes into --out (default results/primer_trimmed_2026-10-06/source_checks_<today>):
#   caproate_trace.csv          one row per study: source cell, unit, value, value in mM C, matrix value
#   caproate_trace_output.txt   everything printed, with the inputs' SHA-256
# Nothing is overwritten: the script stops if either file is already there.

import argparse
import csv
import datetime
import hashlib
import sys
from pathlib import Path

import openpyxl
from openpyxl.utils import column_index_from_string

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "python"))
from build_results_reference import RUN_OF_RECORD, sha256_lf  # noqa: E402

MATRIX_NAME = "BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"
M_CAPROATE, NC_CAPROATE = 116.16, 6   # molar mass and carbon number used for every study
OUT_CSV, OUT_TXT = "caproate_trace.csv", "caproate_trace_output.txt"

# Paper_ID, Sample_ID, file, sheet, caproate column, header row, data row, source unit
TRACES = [
    ("Duber_2020", "ERR3200168", "2020_Evaluation.....xlsx", "UASB II mmol C", "O", 4, 82, "mM C"),
    ("Duber_2022", "SRR18079352", "2022_Lactate and ethanol....xlsx", "LA", "AK", 2, 62, "g/L"),
    ("Duber_2024", "SRR28047211", "2024_Single vs two stage....xlsx", "R1 mmol C", "O", 4, 159, "mM C"),
    ("Brodowski_2022", "SRR18962132", "Brodowski STOTEN2022 Influence of lactate to acetate ratio.xlsx", "B2", "J", 2, 68,
     "mM C"),
    ("Brodowski_2022_EXTERNAL_ACETATE", "SRR13989643", "Brodowski STOTEN2021 Effect of external acetate.xlsx", "B1", "J",
     2, 91, "mM C"),
    ("Brodowski_2025", "SRR28044222", "Brodowski CEJ 2025.xlsx", "BLE", "J", 2, 57, "mM C"),
    ("Hanna_2025", "sample_B4", "metadata_Expanding the product_Prusak.xlsx", "SUCC_META", "Q", 1, 20, "g/L"),
]
# The column each source row's day is read from
DAY_COL = {"Duber_2020": "E", "Duber_2022": "D", "Duber_2024": "E", "Brodowski_2022": "A",
           "Brodowski_2022_EXTERNAL_ACETATE": "A", "Brodowski_2025": "A", "Hanna_2025": "E"}

ap = argparse.ArgumentParser(description="Trace one caproate value per study from the authors' files to the matrix.")
ap.add_argument("authors_dir", type=Path, help="folder holding the seven authors' workbooks")
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


say("Inputs (SHA-256; the matrix with line endings read as LF, the workbooks as saved):")
say(f"  {shown(FINAL)}  {sha256_lf(FINAL)}")
for fname in sorted({t[2] for t in TRACES}):
    say(f"  {fname}  {hashlib.sha256((args.authors_dir / fname).read_bytes()).hexdigest()}")
say()

with open(FINAL, encoding="utf-8-sig", newline="") as f:
    final = {(r["Paper_ID"], r["Sample_ID"]): r for r in csv.DictReader(f)}
out = []
for pid, sid, fname, sheet, col, hrow, drow, unit in TRACES:
    ws = openpyxl.load_workbook(args.authors_dir / fname, data_only=True, read_only=True)[sheet]
    rows = list(ws.iter_rows(values_only=True))
    j = column_index_from_string(col) - 1
    header = rows[hrow - 1][j]
    raw = float(rows[drow - 1][j])
    day = rows[drow - 1][column_index_from_string(DAY_COL[pid]) - 1]
    mmc = raw * 1000 / (M_CAPROATE / NC_CAPROATE) if unit == "g/L" else raw
    f = final[(pid, sid)]
    ours = float(f["Caproate"])
    out.append({
        "Paper_ID": pid, "Sample_ID": sid, "reactor_or_trial": f["BIOREACTOR"], "day_in_matrix": f["DAY"],
        "source_file": fname, "sheet": sheet, "column": col, "column_header": header,
        "excel_row": drow, "day_in_source": day, "source_unit": unit, "value_in_source": raw,
        "conversion": "mM C = g/L x 1000 / (116.16/6)" if unit == "g/L" else "none (file in mM C)",
        "source_in_mM_C": round(mmc, 3), "final_matrix_mM_C": ours,
        "ratio_matrix_to_source": round(ours / mmc, 4) if mmc else "",
        "agrees_within_0.51_mM_C": abs(ours - mmc) <= 0.51,
    })
OUT.mkdir(parents=True, exist_ok=True)
with open(OUT / OUT_CSV, "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), lineterminator="\n")
    w.writeheader()
    w.writerows(out)
for r in out:
    say(f"{r['Paper_ID']} {r['Sample_ID']}: {r['source_file']}, sheet {r['sheet']}, column {r['column']} "
        f"({r['column_header']!r}), row {r['excel_row']}, day {r['day_in_source']}: {r['value_in_source']} {r['source_unit']}"
        f" = {r['source_in_mM_C']} mM C; matrix {r['final_matrix_mM_C']}; agrees: {r['agrees_within_0.51_mM_C']}")
say(f"\nAgree within 0.51 mM C: {sum(r['agrees_within_0.51_mM_C'] for r in out)} of {len(out)}")
say(f"Wrote {shown(OUT / OUT_CSV)}")
(OUT / OUT_TXT).write_text("\n".join(LOG) + "\n", encoding="utf-8")
print(f"Wrote {shown(OUT / OUT_TXT)}")
