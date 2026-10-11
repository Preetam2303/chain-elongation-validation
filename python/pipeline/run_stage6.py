# run_stage6.py
# Stage 6 of the primer-trimmed rerun: run the downstream analyses (02 to 10)
# on the corrected run's matrices, without touching data/ or results/.
#
# Each script reads its matrix from an environment variable when it is set,
# and from data/ otherwise. This runner sets them to the files in the run
# folder's corrected_<date> folder (made by 04_metadata_merge.R, 00 and 00b):
#   BIOTWIN_MATRIX       BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv  (00b; 02-06, 08-10)
#   BIOTWIN_PRUNED       BIOTWIN_PRUNED_ML_MATRIX.csv              (00; 07a, 07b, 10)
#   BIOTWIN_GENUS        BIOTWIN_GENUS_ML_MATRIX.csv               (00; 07c, 07d, 09)
#   BIOTWIN_RESULTS_DIR  the new dated output folder               (02-06, 09, 10 write tables)
#
# Output: a new folder stage6_corrected_<date> inside the corrected_<date>
# folder, with one <script>.log per script (everything it printed), the tables
# the scripts write, and versions.txt. It refuses to run if that folder already
# exists.
#
# Usage, from the repository root:
#   python python/pipeline/run_stage6.py            all scripts, in order
#   python python/pipeline/run_stage6.py 02 04 09   only the scripts named
# 03 (the 999-shuffle permutation test) takes the longest.
#
# The published numbers (and the 2026-10-10 published-settings check run,
# with its --published environment check) come from the code at commit
# 844f478; the settings changed after that.

import datetime
import glob
import importlib.metadata
import os
import platform
import subprocess
import sys

# ============================================================================
# EDIT THIS: the folder that holds your primer_trimmed_* run folder.
# ============================================================================
GRAND_MERGE_DIR = "C:/Users/IISiIS-ZWW-233/Documents/BIOTWIN_Grand_Merge"

SCRIPTS = [
    "02_loso_by_study.py",
    "03_permutation_test_999.py",
    "04_vessel_grouped_deployability.py",
    "05_intrastudy_transferability_table5.py",
    "06_cross_platform_transfer.py",
    "07a_algorithm_comparison_lightgbm_rf.py",
    "07b_algorithm_comparison_lightgbm_standard.py",
    "07c_figure3_intermediate_unscaled_xgboost.py",
    "07d_figure3_stage1_random_forest.py",
    "08_shared_genus_transfer_control.py",
    "09_guild_membership_shap_ranks.py",
    "10_table3a_naive_splits.py",
]
PIPELINE_DIR = os.path.dirname(os.path.abspath(__file__))

wanted = sys.argv[1:]
scripts = [s for s in SCRIPTS if not wanted or s.split("_")[0] in wanted]
unknown = [w for w in wanted if w not in {s.split("_")[0] for s in SCRIPTS}]
if unknown or not scripts:
    raise SystemExit(f"Unknown script number(s): {unknown}. Choose from: {[s.split('_')[0] for s in SCRIPTS]}")

run_dirs = glob.glob(os.path.join(GRAND_MERGE_DIR, "primer_trimmed_*"))
if len(run_dirs) != 1:
    raise SystemExit(f"Expected exactly one primer_trimmed_* folder in {GRAND_MERGE_DIR}, found {len(run_dirs)}")
RUN_DIR = run_dirs[0]
corrected_dirs = glob.glob(os.path.join(RUN_DIR, "corrected_????-??-??"))
if len(corrected_dirs) != 1:
    raise SystemExit(f"Expected exactly one corrected_<date> folder in {RUN_DIR}, found {len(corrected_dirs)}")
CORRECTED_DIR = corrected_dirs[0]
inputs = {
    "BIOTWIN_MATRIX": os.path.join(CORRECTED_DIR, "BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"),
    "BIOTWIN_PRUNED": os.path.join(CORRECTED_DIR, "BIOTWIN_PRUNED_ML_MATRIX.csv"),
    "BIOTWIN_GENUS": os.path.join(CORRECTED_DIR, "BIOTWIN_GENUS_ML_MATRIX.csv"),
}
inputs = {k: os.path.abspath(v) for k, v in inputs.items()}
for name, path in inputs.items():
    if not os.path.exists(path):
        raise SystemExit(f"Missing input {path} ({name}). Run 04_metadata_merge.R, 00 and 00b first.")

OUT_DIR = os.path.join(CORRECTED_DIR, f"stage6_corrected_{datetime.date.today():%Y-%m-%d}")
if os.path.exists(OUT_DIR):
    raise SystemExit(f"{OUT_DIR} already exists. Rename or move it before rerunning.")
os.makedirs(OUT_DIR)
print(f"Corrected-run folder: {CORRECTED_DIR}")
print(f"Writing to: {OUT_DIR}")

env = dict(os.environ, **inputs, BIOTWIN_RESULTS_DIR=OUT_DIR,
           PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8")

with open(os.path.join(OUT_DIR, "versions.txt"), "w", encoding="utf-8") as f:
    f.write(f"Date: {datetime.datetime.now():%Y-%m-%d %H:%M}\n")
    f.write(f"Python {platform.python_version()} on {platform.platform()}\n")
    for pkg in ("pandas", "numpy", "scipy", "scikit-learn", "xgboost", "lightgbm", "shap"):
        try:
            f.write(f"{pkg}=={importlib.metadata.version(pkg)}\n")
        except importlib.metadata.PackageNotFoundError:
            f.write(f"{pkg}: not installed\n")
    for name, path in inputs.items():
        f.write(f"{name}={path}\n")
    f.write(f"Scripts: {', '.join(scripts)}\n")

failed = []
for script in scripts:
    log_path = os.path.join(OUT_DIR, script.replace(".py", ".log"))
    print(f"\n========== {script} ==========")
    with open(log_path, "w", encoding="utf-8") as log:
        proc = subprocess.Popen([sys.executable, script], cwd=PIPELINE_DIR, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, encoding="utf-8", errors="replace")
        for line in proc.stdout:
            print(line, end="")
            log.write(line)
        proc.wait()
    if proc.returncode != 0:
        failed.append(script)
        print(f"!! {script} failed (exit code {proc.returncode}); see {log_path}")

print("\n========== Stage 6 done ==========")
print(f"Logs and outputs: {OUT_DIR}")
if failed:
    print(f"Failed: {', '.join(failed)}")
    sys.exit(1)
