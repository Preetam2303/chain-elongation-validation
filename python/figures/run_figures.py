# run_figures.py
# Stage 8 of the primer-trimmed rerun: draw every figure from the corrected
# run's results files, without touching data/, results/ or figures_output/.
#
# The figures read the tables that Stage 6 wrote (02 to 10), 09's log for its
# fold count, and the corrected genus matrix (00b) for Figures 4, 5, 6, S4 and
# S5. Each Stage 6 script's run of record is the one stage6_corrected_* folder
# inside the corrected_<date> folder that holds its log (<script>_*.log), as in
# build_results_reference.py. After one full run_stage6.py run that is the same
# folder for every script; in the run of 2026-10-11, where some scripts were
# rerun later, it is a later folder for those. If a log is in more than one
# folder, move the one that is not the run of record out of corrected_<date>.
#
# Output: a new folder figures_corrected_<date> inside the corrected_<date>
# folder, with each figure as PNG and PDF, one <script>.log per figure script
# (everything it printed, including every value it plots and each input file
# with its SHA-256), and versions.txt. When only some figures are drawn, the
# folder name ends with their names (e.g. figures_corrected_<date>_2_S4). It
# refuses to run if that folder already exists.
#
# Usage, from the repository root:
#   python python/figures/run_figures.py            every figure, in order
#   python python/figures/run_figures.py 2 6 S4     only the figures named

import datetime
import glob
import importlib.metadata
import os
import platform
import subprocess
import sys

# ============================================================================
# EDIT THIS: the folder that holds your primer_trimmed_* run folder (the same
# folder as in python/pipeline/run_stage6.py).
# ============================================================================
GRAND_MERGE_DIR = "C:/Users/IISiIS-ZWW-233/Documents/BIOTWIN_Grand_Merge"

# The Stage 6 scripts whose tables or logs the figures read
STAGE6_SCRIPTS = ["02", "03", "04", "05", "06", "07c", "07d", "09", "10"]

FIGURES = {
    "1": "figure1_r2_comparison.py",
    "2": "figure2_permutation_histogram.py",
    "3": "figure3_duber_progression.py",
    "4": "figure4_cohort_comparison.py",
    "5": "figure5_shap_genus.py",
    "6": "figure6_attribution_vs_spread.py",
    "S1": "figure_s1_substrate_ladder.py",
    "S3": "figure_s3_intrastudy_transfer.py",
    "S4": "figure_s4_duber_scatter.py",
    "S5": "figure_s5_combined_shap.py",
}
FIGURES_DIR = os.path.dirname(os.path.abspath(__file__))

wanted = sys.argv[1:]
unknown = [w for w in wanted if w not in FIGURES]
if unknown:
    raise SystemExit(f"Unknown figure(s): {unknown}. Choose from: {list(FIGURES)}")
figures = [f for f in FIGURES if not wanted or f in wanted]

run_dirs = glob.glob(os.path.join(GRAND_MERGE_DIR, "primer_trimmed_*"))
if len(run_dirs) != 1:
    raise SystemExit(f"Expected exactly one primer_trimmed_* folder in {GRAND_MERGE_DIR}, found {len(run_dirs)}")
corrected_dirs = glob.glob(os.path.join(run_dirs[0], "corrected_????-??-??"))
if len(corrected_dirs) != 1:
    raise SystemExit(f"Expected exactly one corrected_<date> folder in {run_dirs[0]}, found {len(corrected_dirs)}")
CORRECTED_DIR = os.path.abspath(corrected_dirs[0])

inputs = {"BIOTWIN_MATRIX": os.path.join(CORRECTED_DIR, "BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")}
if not os.path.exists(inputs["BIOTWIN_MATRIX"]):
    raise SystemExit(f"Missing {inputs['BIOTWIN_MATRIX']}. Run 04_metadata_merge.R, 00 and 00b first.")
for script in STAGE6_SCRIPTS:
    found = [d for d in sorted(glob.glob(os.path.join(CORRECTED_DIR, "stage6_corrected_*")))
             if glob.glob(os.path.join(d, f"{script}_*.log"))]
    if len(found) != 1:
        raise SystemExit(f"Expected one stage6_corrected_* folder in {CORRECTED_DIR} with {script}'s log "
                         f"({script}_*.log), found {len(found)}: {[os.path.basename(d) for d in found]}. "
                         + ("Run python/pipeline/run_stage6.py first." if not found else
                            "Move the folder that is not its run of record out of the corrected_<date> folder."))
    inputs[f"BIOTWIN_STAGE6_{script}"] = found[0]

only = "".join("_" + f for f in figures) if wanted else ""
OUT_DIR = os.path.join(CORRECTED_DIR, f"figures_corrected_{datetime.date.today():%Y-%m-%d}{only}")
if os.path.exists(OUT_DIR):
    raise SystemExit(f"{OUT_DIR} already exists. Rename or move it before rerunning.")
os.makedirs(OUT_DIR)
print(f"Corrected-run folder: {CORRECTED_DIR}")
for script in STAGE6_SCRIPTS:
    print(f"Stage 6 {script:<3} from: {os.path.basename(inputs[f'BIOTWIN_STAGE6_{script}'])}")
print(f"Writing to: {OUT_DIR}")

env = dict(os.environ, **inputs, BIOTWIN_FIGURES_DIR=OUT_DIR, MPLBACKEND="Agg",
           PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8")

with open(os.path.join(OUT_DIR, "versions.txt"), "w", encoding="utf-8") as f:
    f.write(f"Date: {datetime.datetime.now():%Y-%m-%d %H:%M}\n")
    f.write(f"Python {platform.python_version()} on {platform.platform()}\n")
    for pkg in ("pandas", "numpy", "scipy", "scikit-learn", "xgboost", "shap", "matplotlib", "seaborn"):
        try:
            f.write(f"{pkg}=={importlib.metadata.version(pkg)}\n")
        except importlib.metadata.PackageNotFoundError:
            f.write(f"{pkg}: not installed\n")
    for name, path in inputs.items():
        f.write(f"{name}={path}\n")
    f.write(f"Figures: {', '.join(figures)}\n")

failed = []
for fig in figures:
    script = FIGURES[fig]
    log_path = os.path.join(OUT_DIR, script.replace(".py", ".log"))
    print(f"\n========== Figure {fig}: {script} ==========")
    with open(log_path, "w", encoding="utf-8") as log:
        proc = subprocess.Popen([sys.executable, script], cwd=FIGURES_DIR, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, encoding="utf-8", errors="replace")
        for line in proc.stdout:
            print(line, end="")
            log.write(line)
        proc.wait()
    if proc.returncode != 0:
        failed.append(script)
        print(f"!! {script} failed (exit code {proc.returncode}); see {log_path}")

print("\n========== Figures done ==========")
print(f"Figures and logs: {OUT_DIR}")
if failed:
    print(f"Failed: {', '.join(failed)}")
    sys.exit(1)
