# 03b_reactor_aware_permutation.py
# Supplementary sensitivity for the permutation controls of 03 (Section 3.4,
# Figure 2), added 2026-10-19 after the 15 Oct critic review (comment 1). 03's
# within-study shuffle moves caproate between the reactors of a study and
# between sampling days, so its null treats every sampling unit of a study as
# exchangeable. This script refits the same headline LOSO (full matrix, core
# predictor set, the same folds and settings as 03) under two shuffles that
# keep the reactor structure. 03 and its outputs are unchanged; both shuffles
# are reported whatever they show (decision 2026-10-19).
#   1. Within reactor: caproate is permuted only among the sampling units of
#      the same reactor or bottle label (Paper_ID + BIOREACTOR), so each
#      reactor keeps its own values.
#   2. Cyclic shift within reactor: each reactor's caproate series, in order of
#      sampling day, is shifted by a random lag of 1 to n - 1 units (n = the
#      reactor's sampling units), so each reactor keeps its values and their
#      autocorrelation in time; only the alignment of caproate with that day's
#      predictors is broken.
# Both move caproate by sampling unit, as 03 does: the same-day replicates of
# Duber 2024 and the Nanopore study share one value and move together.
# Reactors with a single sampling unit (Brodowski 2022 B1, R2, R3 and R4) keep
# their value. Each Duber 2022 bottle label (R1-R9) holds two bottles sampled
# on day 22; both shuffles can only swap the two, and the cyclic shift always
# does (a lag of 1 of 2).
#
# Shuffle i uses np.random.RandomState(i), drawing the reactors in
# alphabetical order of their label. Sampling days are read from DAY with the
# Nanopore study's leading "d" removed; units of one reactor on the same day
# keep the order of their sampling-unit id. p = (k + 1) / (N + 1), as in 03:
# k counts the shuffles whose mean R2 reaches the observed one (rounded to 3
# decimals), or whose mean RMSE is as low as the observed one or lower. The
# same p is also given without the Duber 2022 fold (the mean over the other
# six held-out studies, from the same shuffles), because most of the gain over
# the training-mean benchmark comes from that study.
#
# Output (in BIOTWIN_RESULTS_DIR, default ../../results):
#   03b_permutation_null_within_reactor_<N>.csv   mean LOSO R2 and RMSE per shuffle,
#                                                 then R2 and RMSE per held-out study
#   03b_permutation_null_cyclic_shift_<N>.csv     the same for the cyclic shift
#   03b_reactor_aware_permutation_tests.csv       one row per shuffle: observed, null, k and p
#   03b_reactor_aware_permutation_tests_by_study.csv   the same per held-out study
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
LEFT_OUT_STUDY = 'Duber_2022'
NULL_PATHS = {
    'within reactor': os.path.join(OUTPUT_DIR, f'03b_permutation_null_within_reactor_{N_ITERATIONS}.csv'),
    'cyclic shift within reactor': os.path.join(OUTPUT_DIR, f'03b_permutation_null_cyclic_shift_{N_ITERATIONS}.csv'),
}
SUMMARY_PATH = os.path.join(OUTPUT_DIR, '03b_reactor_aware_permutation_tests.csv')
BY_STUDY_PATH = os.path.join(OUTPUT_DIR, '03b_reactor_aware_permutation_tests_by_study.csv')
for path in list(NULL_PATHS.values()) + [SUMMARY_PATH, BY_STUDY_PATH]:
    if os.path.exists(path):  # checked before the long run, not after it
        raise SystemExit(f"{path} already exists. Set BIOTWIN_RESULTS_DIR to a new folder.")

df, genus_cols = load_matrix(MATRIX_PATH)
features = sorted(set(present(OPS + CORE_CHEM + CAT_COLS, df) + genus_cols))
X = df[features]
papers = df['Paper_ID']
y_true = df[TARGET]

# The folds of 03, with their scaled predictors computed once; only y changes per shuffle.
folds, held_out = [], []
for train_idx, test_idx in LeaveOneGroupOut().split(X, y_true, groups=papers):
    if len(test_idx) > 1:
        X_train, X_test = scale(X.iloc[train_idx], X.iloc[test_idx])
        folds.append((train_idx, test_idx, X_train, X_test))
        held_out.append(papers.iloc[test_idx[0]])
if LEFT_OUT_STUDY not in held_out:
    raise SystemExit(f"{LEFT_OUT_STUDY} is not one of the held-out studies {held_out}.")
kept = np.array([h != LEFT_OUT_STUDY for h in held_out])


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
unit_reactor = (df['Paper_ID'] + '|' + df['BIOREACTOR']).groupby(units).first().reindex(unit_order)
day_text = df['DAY'].astype(str).str.strip().str.removeprefix('d')
readable = day_text.str.fullmatch(r'\d+(\.\d+)?')
if not readable.all():
    raise SystemExit(f"Unreadable sampling days: {sorted(df.loc[~readable, 'DAY'].astype(str).unique())}")
unit_day = day_text.astype(float).groupby(units).first().reindex(unit_order)
REACTORS = sorted(unit_reactor.unique())
# Each reactor's sampling units: in id order (within reactor) and in day order (cyclic shift).
REACTOR_UNITS = {r: unit_order[unit_reactor.values == r] for r in REACTORS}
REACTOR_UNITS_BY_DAY = {r: idx[np.argsort(unit_day.loc[idx].values, kind='mergesort')]
                        for r, idx in REACTOR_UNITS.items()}


def to_rows(shuffled_units):
    return pd.Series(units.map(shuffled_units).values, index=df.index)


def shuffled_within_reactor(seed):
    rng = np.random.RandomState(seed)
    shuffled = unit_values.copy()
    for reactor in REACTORS:
        idx = REACTOR_UNITS[reactor]
        shuffled.loc[idx] = rng.permutation(unit_values.loc[idx].values)
    return to_rows(shuffled)


def shifted_within_reactor(seed):
    rng = np.random.RandomState(seed)
    shuffled = unit_values.copy()
    for reactor in REACTORS:
        idx = REACTOR_UNITS_BY_DAY[reactor]
        if len(idx) > 1:
            shuffled.loc[idx] = np.roll(unit_values.loc[idx].values, rng.randint(1, len(idx)))
    return to_rows(shuffled)


observed_folds_r2, observed_folds_rmse = loso_r2_rmse(y_true)
observed_r2, observed_rmse = observed_folds_r2.mean(), observed_folds_rmse.mean()
TRUE_OBSERVED_R2 = round(observed_r2, 3)
observed_r2_kept = round(observed_folds_r2[kept].mean(), 3)
observed_rmse_kept = observed_folds_rmse[kept].mean()
benchmark_rmse = np.mean([training_mean_r2_rmse(y_true.iloc[tr], y_true.iloc[te])[1] for tr, te, _, _ in folds])

print(f"\n--- {N_ITERATIONS}-SHUFFLE REACTOR-AWARE PERMUTATION TESTS (headline: full-matrix LOSO, core set) ---")
print(f"Rows: {len(df)} | sampling units: {len(unit_values)} in {len(REACTORS)} reactors "
      f"(same-day replicates move together) | predictors: {len(features)}")
print(f"Observed LOSO R2 (true labels): {TRUE_OBSERVED_R2} | observed mean RMSE: {observed_rmse:.2f} mM C "
      f"(training-mean benchmark: {benchmark_rmse:.1f})")
print(f"Without the {LEFT_OUT_STUDY} fold: observed R2 {observed_r2_kept}, RMSE {observed_rmse_kept:.2f} mM C")
print("Sampling units per reactor, in day order (single-unit reactors keep their value):")
for reactor in REACTORS:
    idx = REACTOR_UNITS_BY_DAY[reactor]
    print(f"  {reactor:40s} {len(idx):2d} unit(s), days " + ", ".join(f"{d:g}" for d in unit_day.loc[idx]))

summary, by_study = [], []
for scheme, shuffle in [('within reactor', shuffled_within_reactor),
                        ('cyclic shift within reactor', shifted_within_reactor)]:
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
    per = pd.DataFrame(per_study)
    r2s_kept = per[[f'r2_{h}' for h in np.array(held_out)[kept]]].mean(axis=1)
    rmses_kept = per[[f'rmse_{h}' for h in np.array(held_out)[kept]]].mean(axis=1)
    k_r2 = sum(r2 >= TRUE_OBSERVED_R2 for r2 in r2s)
    k_rmse = sum(rmse <= observed_rmse for rmse in rmses)
    k_r2_kept = int((r2s_kept >= observed_r2_kept).sum())
    k_rmse_kept = int((rmses_kept <= observed_rmse_kept).sum())
    summary.append({'shuffle': scheme, 'n_shuffles': N_ITERATIONS, 'observed_r2': TRUE_OBSERVED_R2,
                    'null_r2_mean': np.mean(r2s), 'null_r2_sd': np.std(r2s), 'null_r2_median': np.median(r2s),
                    'k_r2': k_r2, 'p_r2': (k_r2 + 1) / (N_ITERATIONS + 1),
                    'observed_rmse': observed_rmse, 'null_rmse_mean': np.mean(rmses),
                    'null_rmse_sd': np.std(rmses), 'k_rmse': k_rmse, 'p_rmse': (k_rmse + 1) / (N_ITERATIONS + 1),
                    'training_mean_benchmark_rmse': benchmark_rmse,
                    'left_out_study': LEFT_OUT_STUDY, 'observed_r2_without': observed_r2_kept,
                    'k_r2_without': k_r2_kept, 'p_r2_without': (k_r2_kept + 1) / (N_ITERATIONS + 1),
                    'observed_rmse_without': observed_rmse_kept,
                    'k_rmse_without': k_rmse_kept, 'p_rmse_without': (k_rmse_kept + 1) / (N_ITERATIONS + 1)})
    for i, h in enumerate(held_out):
        study_r2s = per[f'r2_{h}']
        study_rmses = per[f'rmse_{h}']
        obs_r2 = round(observed_folds_r2[i], 3)
        k_r2_h = int((study_r2s >= obs_r2).sum())
        k_rmse_h = int((study_rmses <= observed_folds_rmse[i]).sum())
        by_study.append({'shuffle': scheme, 'held_out_study': h, 'n_shuffles': N_ITERATIONS,
                         'observed_r2': obs_r2, 'null_r2_mean': study_r2s.mean(),
                         'null_r2_median': study_r2s.median(),
                         'k_r2': k_r2_h, 'p_r2': (k_r2_h + 1) / (N_ITERATIONS + 1),
                         'observed_rmse': observed_folds_rmse[i], 'null_rmse_mean': study_rmses.mean(),
                         'k_rmse': k_rmse_h, 'p_rmse': (k_rmse_h + 1) / (N_ITERATIONS + 1)})
    print(f"  Saved full results to: {NULL_PATHS[scheme]}")

print("-" * 80)
print(f"NULL DISTRIBUTIONS (N={N_ITERATIONS}), observed R2 {TRUE_OBSERVED_R2:.3f}, RMSE {observed_rmse:.1f} mM C:")
for s in summary:
    print(f"  Shuffled {s['shuffle']:28s} R2 {s['null_r2_mean']:7.3f} (± {s['null_r2_sd']:.3f}), "
          f"median {s['null_r2_median']:7.3f} | k = {s['k_r2']:3d}, p = {s['p_r2']:.4f} || "
          f"RMSE {s['null_rmse_mean']:6.1f} (± {s['null_rmse_sd']:.1f}) | k = {s['k_rmse']:3d}, p = {s['p_rmse']:.4f}")
    print(f"    without the {LEFT_OUT_STUDY} fold: p(R2) {s['p_r2_without']:.4f}, p(RMSE) {s['p_rmse_without']:.4f}")
print("  k = shuffles that do at least as well as the true labels (R2 as high, or RMSE as low)")
print("Per held-out study (observed R2 and RMSE; p of each shuffle scheme):")
for i, h in enumerate(held_out):
    rows = [b for b in by_study if b['held_out_study'] == h]
    print(f"  {h:34s} R2 {rows[0]['observed_r2']:7.3f}, RMSE {rows[0]['observed_rmse']:6.1f} | " +
          " | ".join(f"{b['shuffle']}: p(R2) {b['p_r2']:.4f}, p(RMSE) {b['p_rmse']:.4f}" for b in rows))
pd.DataFrame(summary).to_csv(SUMMARY_PATH, index=False)
pd.DataFrame(by_study).to_csv(BY_STUDY_PATH, index=False)
print(f"  Saved: {SUMMARY_PATH}")
print(f"  Saved: {BY_STUDY_PATH}")
print("-" * 80)
