# 03_permutation_test_999.py
# Target-permutation controls for the headline (Section 3.4, Figure 2): the
# full-matrix LOSO with exactly the core predictor set, the same model as
# Table 3b's LOSO row and Table 4's full-matrix tier 3 (02_loso_by_study.py).
#
# Two shuffles of caproate, reported side by side (decision 2026-10-11):
#   1. Across studies: caproate is shuffled over the sampling units of all
#      studies. This mixes values between studies, which erases the
#      between-study differences that make the real task hard.
#   2. Within each study: caproate is shuffled among the sampling units of the
#      same study only, so every study keeps its own values and the
#      between-study differences stay. This asks whether the predictors (the
#      genera and the chemistry together) carry information that transfers
#      beyond study-level patterns. It does not keep the reactor a value came
#      from or the order of sampling days; 03b_reactor_aware_permutation.py
#      adds two shuffles that do (exploratory sensitivity analyses).
#
# Both shuffle by sampling unit, not by row (decision 2026-10-08): the three
# same-day replicate samples of Duber 2024 and of the Nanopore study share one
# caproate value, so they move together; every other row is its own unit.
# Shuffle i uses random_state = i (across studies, unchanged) or
# np.random.RandomState(i) drawing the studies in alphabetical order (within
# each study). The observed value is computed in the same run, with the true
# labels and the same features and folds.
#
# Each shuffle records the mean LOSO R2 (the headline statistic) and the mean
# LOSO RMSE, and both per held-out study. p = (k + 1) / (N + 1), where k is
# the number of shuffles whose R2 is as high as the observed one or higher, or,
# for RMSE, whose RMSE is as low as the observed one or lower. Every comparison and every
# saved statistic is at full precision; values are rounded only when printed
# (since 2026-10-19: before, R2 was compared with the observed value rounded
# to 3 decimals, which moved two per-study k by R2 in the run of record of
# 2026-10-11; build_results_reference.py recomputes them from the saved
# shuffles).
#
# Self-check: the within-study shuffle keeps each study's set of caproate
# values when every sampling unit of the study has the same number of rows
# (checked here), so for each held-out study the shuffles that reach its
# observed R2 are exactly those that reach its observed RMSE, and p by R2 must
# equal p by RMSE. The script prints any mismatch. The across-study shuffle
# changes each study's values, so it has no such check.
#
# Output (in BIOTWIN_RESULTS_DIR, default ../../results):
#   permutation_null_distribution_<N>.csv   across-study shuffles: mean LOSO R2 and RMSE,
#                                           then R2 and RMSE per held-out study
#   permutation_null_within_study_<N>.csv   within-study shuffles: the same
#   03_permutation_tests.csv                one row per shuffle scheme: observed, null, k and p
#   03_permutation_tests_by_study.csv       the same per held-out study
# Set BIOTWIN_N_PERMUTATIONS (default 999) to e.g. 20 for a quick check.

import os
import time
import warnings

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import LeaveOneGroupOut

from pipeline_settings import (CAT_COLS, CORE_CHEM, OPS, TARGET, XGB_PARAMS, load_matrix, present, r2_rmse,
                               sampling_units, scale, training_mean_r2_rmse)

warnings.filterwarnings('ignore')

MATRIX_PATH = os.environ.get("BIOTWIN_MATRIX", "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")
OUTPUT_DIR = os.environ.get("BIOTWIN_RESULTS_DIR", "../../results")
os.makedirs(OUTPUT_DIR, exist_ok=True)
N_ITERATIONS = int(os.environ.get("BIOTWIN_N_PERMUTATIONS", 999))
NULL_PATHS = {
    'across studies': os.path.join(OUTPUT_DIR, f'permutation_null_distribution_{N_ITERATIONS}.csv'),
    'within each study': os.path.join(OUTPUT_DIR, f'permutation_null_within_study_{N_ITERATIONS}.csv'),
}
SUMMARY_PATH = os.path.join(OUTPUT_DIR, '03_permutation_tests.csv')
BY_STUDY_PATH = os.path.join(OUTPUT_DIR, '03_permutation_tests_by_study.csv')
for path in list(NULL_PATHS.values()) + [SUMMARY_PATH, BY_STUDY_PATH]:
    if os.path.exists(path):  # checked before the long run, not after it
        raise SystemExit(f"{path} already exists. Set BIOTWIN_RESULTS_DIR to a new folder.")

df, genus_cols = load_matrix(MATRIX_PATH)
features = sorted(set(present(OPS + CORE_CHEM + CAT_COLS, df) + genus_cols))
X = df[features]
papers = df['Paper_ID']
y_true = df[TARGET]

# Precompute the folds and their scaled predictors once; only y changes per shuffle.
folds, held_out = [], []
for train_idx, test_idx in LeaveOneGroupOut().split(X, y_true, groups=papers):
    if len(test_idx) > 1:
        X_train, X_test = scale(X.iloc[train_idx], X.iloc[test_idx])
        folds.append((train_idx, test_idx, X_train, X_test))
        held_out.append(papers.iloc[test_idx[0]])


def loso_r2_rmse(y):
    """R2 and RMSE of each held-out study, in the order of held_out."""
    fold_r2, fold_rmse = [], []
    for train_idx, test_idx, X_train, X_test in folds:
        model = xgb.XGBRegressor(**XGB_PARAMS).fit(X_train, y.iloc[train_idx])
        r2, rmse = r2_rmse(y.iloc[test_idx], model.predict(X_test))
        fold_r2.append(r2); fold_rmse.append(rmse)
    return np.array(fold_r2), np.array(fold_rmse)


units = sampling_units(df)
unit_values = y_true.groupby(units).first()   # one caproate value per sampling unit
unit_order = unit_values.index
unit_study = papers.groupby(units).first().reindex(unit_order)
STUDIES = sorted(unit_study.unique())
# The self-check's condition: every sampling unit of a study has the same number of rows.
rows_per_unit = units.value_counts().groupby(unit_study).nunique()
UNEVEN_STUDIES = sorted(rows_per_unit.index[rows_per_unit > 1])


def to_rows(shuffled_units):
    return pd.Series(units.map(shuffled_units).values, index=df.index)


def shuffled_across(seed):
    return to_rows(pd.Series(unit_values.sample(frac=1.0, random_state=seed).values, index=unit_order))


def shuffled_within(seed):
    rng = np.random.RandomState(seed)
    shuffled = unit_values.copy()
    for study in STUDIES:
        idx = unit_order[unit_study.values == study]
        shuffled.loc[idx] = rng.permutation(unit_values.loc[idx].values)
    return to_rows(shuffled)


observed_folds_r2, observed_folds_rmse = loso_r2_rmse(y_true)
observed_r2, observed_rmse = observed_folds_r2.mean(), observed_folds_rmse.mean()
benchmark_rmse = np.mean([training_mean_r2_rmse(y_true.iloc[tr], y_true.iloc[te])[1] for tr, te, _, _ in folds])

print(f"\n--- {N_ITERATIONS}-SHUFFLE TARGET PERMUTATION TESTS (headline: full-matrix LOSO, core set) ---")
print(f"Rows: {len(df)} | sampling units: {len(unit_values)} in {len(STUDIES)} studies "
      f"(same-day replicates move together) | predictors: {len(features)}")
print(f"Observed LOSO R2 (true labels): {observed_r2:.3f} | observed mean RMSE: {observed_rmse:.1f} mM C "
      f"(training-mean benchmark: {benchmark_rmse:.1f})")

summary, by_study = [], []
for scheme, shuffle in [('across studies', shuffled_across), ('within each study', shuffled_within)]:
    print("-" * 80)
    print(f"Shuffling caproate {scheme}")
    r2s, rmses, per_study = [], [], []
    start_time = time.time()
    for iteration in range(1, N_ITERATIONS + 1):
        fold_r2, fold_rmse = loso_r2_rmse(shuffle(iteration))
        r2, rmse = fold_r2.mean(), fold_rmse.mean()
        r2s.append(r2); rmses.append(rmse)
        per_study.append({**{f'r2_{h}': v for h, v in zip(held_out, fold_r2)},
                          **{f'rmse_{h}': v for h, v in zip(held_out, fold_rmse)}})
        if iteration % 10 == 0 or iteration == 1 or iteration == N_ITERATIONS:
            elapsed = time.time() - start_time
            print(f"Iteration {iteration:03d}/{N_ITERATIONS} | Permuted LOGO R2: {r2:6.3f} | "
                  f"RMSE: {rmse:6.1f} | Elapsed Time: {elapsed:.1f}s")
    pd.concat([pd.DataFrame({'iteration': range(1, N_ITERATIONS + 1), 'permuted_r2': r2s,
                             'permuted_rmse': rmses}), pd.DataFrame(per_study)],
              axis=1).to_csv(NULL_PATHS[scheme], index=False)
    k_r2 = sum(r2 >= observed_r2 for r2 in r2s)
    k_rmse = sum(rmse <= observed_rmse for rmse in rmses)
    summary.append({'shuffle': scheme, 'n_shuffles': N_ITERATIONS, 'observed_r2': observed_r2,
                    'null_r2_mean': np.mean(r2s), 'null_r2_sd': np.std(r2s), 'null_r2_median': np.median(r2s),
                    'k_r2': k_r2, 'p_r2': (k_r2 + 1) / (N_ITERATIONS + 1),
                    'observed_rmse': observed_rmse, 'null_rmse_mean': np.mean(rmses),
                    'null_rmse_sd': np.std(rmses), 'k_rmse': k_rmse, 'p_rmse': (k_rmse + 1) / (N_ITERATIONS + 1),
                    'training_mean_benchmark_rmse': benchmark_rmse})
    for i, h in enumerate(held_out):
        study_r2s = [row[f'r2_{h}'] for row in per_study]
        study_rmses = [row[f'rmse_{h}'] for row in per_study]
        obs_r2 = observed_folds_r2[i]
        k_r2_h = sum(r2 >= obs_r2 for r2 in study_r2s)
        k_rmse_h = sum(rmse <= observed_folds_rmse[i] for rmse in study_rmses)
        by_study.append({'shuffle': scheme, 'held_out_study': h, 'n_shuffles': N_ITERATIONS,
                         'observed_r2': obs_r2, 'null_r2_mean': np.mean(study_r2s),
                         'null_r2_median': np.median(study_r2s),
                         'k_r2': k_r2_h, 'p_r2': (k_r2_h + 1) / (N_ITERATIONS + 1),
                         'observed_rmse': observed_folds_rmse[i], 'null_rmse_mean': np.mean(study_rmses),
                         'k_rmse': k_rmse_h, 'p_rmse': (k_rmse_h + 1) / (N_ITERATIONS + 1)})
    print(f"  Saved full results to: {NULL_PATHS[scheme]}")

print("-" * 80)
print(f"NULL DISTRIBUTIONS (N={N_ITERATIONS}), observed R2 {observed_r2:.3f}, RMSE {observed_rmse:.1f} mM C:")
for s in summary:
    print(f"  Shuffled {s['shuffle']:18s} R2 {s['null_r2_mean']:7.3f} (± {s['null_r2_sd']:.3f}), "
          f"median {s['null_r2_median']:7.3f} | k = {s['k_r2']:3d}, p = {s['p_r2']:.4f} || "
          f"RMSE {s['null_rmse_mean']:6.1f} (± {s['null_rmse_sd']:.1f}) | k = {s['k_rmse']:3d}, p = {s['p_rmse']:.4f}")
print("  k = shuffles that do at least as well as the true labels (R2 as high, or RMSE as low)")
print("Per held-out study (observed R2 and RMSE; p of each shuffle scheme):")
for i, h in enumerate(held_out):
    rows = [b for b in by_study if b['held_out_study'] == h]
    print(f"  {h:34s} R2 {rows[0]['observed_r2']:7.3f}, RMSE {rows[0]['observed_rmse']:6.1f} | " +
          " | ".join(f"{b['shuffle']}: p(R2) {b['p_r2']:.4f}, p(RMSE) {b['p_rmse']:.4f}" for b in rows))
print("Self-check (within each study): per held-out study, p by R2 must equal p by RMSE when every sampling "
      "unit of the study has the same number of rows.")
if UNEVEN_STUDIES:
    print(f"  Studies whose sampling units differ in rows (no check): {', '.join(UNEVEN_STUDIES)}")
checked = [b for b in by_study if b['shuffle'] == 'within each study' and b['held_out_study'] not in UNEVEN_STUDIES]
mismatches = [b for b in checked if b['k_r2'] != b['k_rmse']]
for b in mismatches:
    print(f"  MISMATCH {b['held_out_study']}: k(R2) = {b['k_r2']}, k(RMSE) = {b['k_rmse']}")
print(f"  {len(checked) - len(mismatches)} of {len(checked)} held-out studies agree; {len(mismatches)} mismatch(es).")
pd.DataFrame(summary).to_csv(SUMMARY_PATH, index=False)
pd.DataFrame(by_study).to_csv(BY_STUDY_PATH, index=False)
print(f"  Saved: {SUMMARY_PATH}")
print(f"  Saved: {BY_STUDY_PATH}")
print("-" * 80)
