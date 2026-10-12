# figure_inputs.py
# Where the figure scripts read their inputs and write their figures (Stage 8
# of the primer-trimmed rerun). Every figure is drawn from the corrected run's
# results files, the run of record, never from numbers typed into a script.
# run_figures.py sets three kinds of environment variable:
#   BIOTWIN_MATRIX        the corrected genus matrix (00b), for Figures 4, 5, 6, S4 and S5
#   BIOTWIN_STAGE6_<nn>   the folder that holds Stage 6 script <nn>'s run of record
#   BIOTWIN_FIGURES_DIR   the new dated folder the figures are written to
# Every input is printed with its SHA-256, computed with the line endings read
# as LF, so the same file gives the same value on Windows and Linux.
#
# Importing this module also puts python/pipeline/ on the path, so the figure
# scripts can use the shared settings in pipeline_settings.py.

import glob
import hashlib
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))

# How the manuscript names each study (Table 1 and the text)
STUDY_LABELS = {
    "Brodowski_2022": "Brodowski 2022",
    "Brodowski_2022_EXTERNAL_ACETATE": "Brodowski 2022 [Ext.]",
    "Brodowski_2025": "Brodowski 2025",
    "Duber_2020": "Duber 2020",
    "Duber_2022": "Duber 2022",
    "Duber_2024": "Duber 2025",
    "Hanna_2025": "Prusak 2026",
}


def _env(name):
    value = os.environ.get(name)
    if not value:
        raise SystemExit(f"{name} is not set. Run the figures with python/figures/run_figures.py, which points "
                         "them at the corrected run's results files.")
    return value


def sha256_lf(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read().replace(b"\r\n", b"\n")).hexdigest()


def _report(path):
    print(f"Input: {path}\n       sha256 {sha256_lf(path)}")
    return path


def matrix_path():
    return _report(_env("BIOTWIN_MATRIX"))


def _stage6_log(script):
    folder = _env(f"BIOTWIN_STAGE6_{script}")
    logs = glob.glob(os.path.join(folder, f"{script}_*.log"))
    if len(logs) != 1:
        raise SystemExit(f"{folder} should hold one {script}_*.log (the run of {script}); found {len(logs)}.")
    return folder, logs[0]


def results_file(script, name):
    """A table written by Stage 6 script `script`, from its run of record."""
    folder, _ = _stage6_log(script)
    path = os.path.join(folder, name)
    if not os.path.exists(path):
        raise SystemExit(f"Missing {path}: script {script}'s run of record should have written it.")
    return _report(path)


def log_number(script, pattern):
    """The number(s) that `pattern` captures in script's log, for the few
    results that a script prints but does not save. The pattern must match
    exactly one line."""
    _, path = _stage6_log(script)
    _report(path)
    with open(path, encoding="utf-8") as f:
        text = f.read().replace("\r\n", "\n")
    found = re.findall(pattern, text, flags=re.M)
    if len(found) != 1:
        raise SystemExit(f"Expected one line matching {pattern!r} in {path}, found {len(found)}.")
    return tuple(float(x) for x in found[0]) if isinstance(found[0], tuple) else float(found[0])


def save(fig, stem):
    """Write the figure as PNG and PDF into the dated figures folder, never over an existing file."""
    out_dir = _env("BIOTWIN_FIGURES_DIR")
    for ext in ("png", "pdf"):
        path = os.path.join(out_dir, f"{stem}.{ext}")
        if os.path.exists(path):
            raise SystemExit(f"{path} already exists.")
        fig.savefig(path, dpi=300, bbox_inches="tight")
        print(f"Saved: {path}")


def one_row(table, **match):
    """The single row of `table` whose columns equal the given values."""
    rows = table
    for column, value in match.items():
        rows = rows[rows[column] == value]
    if len(rows) != 1:
        raise SystemExit(f"Expected one row with {match}, found {len(rows)}.")
    return rows.iloc[0]
