# 03b_reactor_aware_permutation.py
# Exploratory sensitivity analyses for the permutation controls of 03 (Section
# 3.4, Figure 2), added 2026-10-19 after the 15 Oct critic review (comment 1).
# 03's within-study shuffle moves caproate between the reactors of a study and
# between sampling days, so its null treats every sampling unit of a study as
# exchangeable. This script refits the same headline LOSO (full matrix, core
# predictor set, the same folds and settings as 03) under two shuffles that
# keep caproate inside its reactor. 03's saved outputs are not rewritten; both
# shuffles are reported whatever they show (decisions of 2026-10-19).
#   1. Within reactor: caproate is permuted among the sampling units of the
#      same reactor or bottle label (Paper_ID + BIOREACTOR), so each reactor
#      keeps its own values.
#   2. Cyclic shift within reactor: each reactor's caproate values, in order of
#      sampling day, are rotated by a random lag drawn from 0 to n - 1 (n = the
#      reactor's sampling units; the last values wrap round to the first days).
#      Each reactor keeps its values in the same cyclic order; only their
#      alignment with the days' predictors changes. The sampling intervals are
#      uneven, so this is not a shift by a fixed time.
# Both keep each reactor's values but break their alignment with the days'
# predictors, so a predictor that only tracks time since start-up also counts
# as information, and for a reactor with a strong trend the cyclic shift joins
# its last and first days in one jump. Neither shuffle separates predictive
# skill from a time trend the reactors share.
# Both move caproate by sampling unit, as 03 does: the same-day replicates of
# Duber 2024 and the Nanopore study share one value and move together.
# Reactors with a single sampling unit (Brodowski 2022 B1, R2, R3 and R4) keep
# their value. Each Duber 2022 label (R1-R9) holds two replicate bottles
# sampled once (day 22), so both shuffles can only exchange those two bottles:
# for Duber 2022 either shuffle is a permutation between replicate bottles,
# not a time shift. A reactor of n sampling units has only n cyclic
# arrangements, so under the cyclic shift Duber 2020's held-out values (one
# reactor, 7 units) take only 7 and Brodowski 2022's (only B2, 4 units, can
# move) only 4. In a shuffle that leaves a held-out study's own values in
# place, its R2 and RMSE change only through the shuffled training values of
# the other studies.
#
# Shuffle i uses np.random.RandomState(i), drawing the reactors in sorted
# order of their Paper_ID|BIOREACTOR label. Sampling days are read from DAY
# with any leading "d" removed (only the Nanopore study writes one). Only
# Duber 2022 has two units of one reactor on the same day; they keep the order
# of their sampling-unit id, which for two units does not change the null.
# p = (k + 1) / (N + 1), as in 03: k counts the shuffles whose R2 is as high as
# the observed one or higher, or whose RMSE is as low as the observed one or
# lower. Every comparison and every saved statistic is at full precision;
# values are rounded only when printed. The test is repeated without the Duber
# 2022 held-out fold: the mean over the other six held-out studies, from the
# same shuffles. Duber 2022's samples stay in the training sets of those six
# folds, with their values shuffled like the others'.
#
# Self-check: both shuffles keep each study's set of caproate values when every
# sampling unit of the study has the same number of rows (checked here), so for
# each held-out study the shuffles that reach its observed R2 are exactly those
# that reach its observed RMSE, and p by R2 must equal p by RMSE. The script
# counts, for each study and shuffle scheme, the shuffles that reach one but
# not the other, and prints any mismatch.
#
# Output (in BIOTWIN_RESULTS_DIR, default ../../results):
#   03b_permutation_null_within_reactor_<N>.csv   mean LOSO R2 and RMSE per shuffle,
#                                                 then R2 and RMSE per held-out study
#   03b_permutation_null_cyclic_shift_<N>.csv     the same for the cyclic shift
#   03b_reactor_aware_permutation_tests.csv       one row per shuffle scheme: observed, null, k and p
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
kept = [h for h in held_out if h != LEFT_OUT_STUDY]


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
# The self-check's condition: every sampling unit of a study has the same number of rows.
rows_per_unit = units.value_counts().groupby(papers.groupby(units).first()).nunique()
UNEVEN_STUDIES = sorted(rows_per_unit.index[rows_per_unit > 1])


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
        shuffled.loc[idx] = np.roll(unit_values.loc[idx].values, rng.randint(0, len(idx)))
    return to_rows(shuffled)


observed_folds_r2, observed_folds_rmse = loso_r2_rmse(y_true)
observed = dict(zip(held_out, zip(observed_folds_r2, observed_folds_rmse)))
observed_r2, observed_rmse = observed_folds_r2.mean(), observed_folds_rmse.mean()
observed_r2_kept = np.mean([observed[h][0] for h in kept])
observed_rmse_kept = np.mean([observed[h][1] for h in kept])
benchmark_rmse = np.mean([training_mean_r2_rmse(y_true.iloc[tr], y_true.iloc[te])[1] for tr, te, _, _ in folds])

print(f"\n--- {N_ITERATIONS}-SHUFFLE REACTOR-AWARE PERMUTATION TESTS (exploratory; headline: full-matrix LOSO, "
      "core set) ---")
print(f"Rows: {len(df)} | sampling units: {len(unit_values)} in {len(REACTORS)} reactors "
      f"(same-day replicates move together) | predictors: {len(features)}")
print(f"Observed LOSO R2 (true labels): {observed_r2:.3f} | observed mean RMSE: {observed_rmse:.2f} mM C "
      f"(training-mean benchmark: {benchmark_rmse:.1f})")
print(f"Without the {LEFT_OUT_STUDY} held-out fold: observed R2 {observed_r2_kept:.3f}, "
      f"RMSE {observed_rmse_kept:.2f} mM C")
print("Sampling units per reactor, in day order (single-unit reactors keep their value):")
for reactor in REACTORS:
    idx = REACTOR_UNITS_BY_DAY[reactor]
    print(f"  {reactor:40s} {len(idx):2d} unit(s), days " + ", ".join(f"{d:g}" for d in unit_day.loc[idx]))

summary, by_study, disagree = [], [], {}
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
    r2s, rmses = np.array(r2s), np.array(rmses)
    r2s_kept = per[[f'r2_{h}' for h in kept]].mean(axis=1)
    rmses_kept = per[[f'rmse_{h}' for h in kept]].mean(axis=1)
    k_r2 = int((r2s >= observed_r2).sum())
    k_rmse = int((rmses <= observed_rmse).sum())
    k_r2_kept = int((r2s_kept >= observed_r2_kept).sum())
    k_rmse_kept = int((rmses_kept <= observed_rmse_kept).sum())
    summary.append({'shuffle': scheme, 'n_shuffles': N_ITERATIONS, 'observed_r2': observed_r2,
                    'null_r2_mean': r2s.mean(), 'null_r2_sd': r2s.std(), 'null_r2_median': np.median(r2s),
                    'k_r2': k_r2, 'p_r2': (k_r2 + 1) / (N_ITERATIONS + 1),
                    'observed_rmse': observed_rmse, 'null_rmse_mean': rmses.mean(),
                    'null_rmse_sd': rmses.std(), 'k_rmse': k_rmse, 'p_rmse': (k_rmse + 1) / (N_ITERATIONS + 1),
                    'training_mean_benchmark_rmse': benchmark_rmse,
                    'left_out_study': LEFT_OUT_STUDY, 'observed_r2_without': observed_r2_kept,
                    'k_r2_without': k_r2_kept, 'p_r2_without': (k_r2_kept + 1) / (N_ITERATIONS + 1),
                    'observed_rmse_without': observed_rmse_kept,
                    'k_rmse_without': k_rmse_kept, 'p_rmse_without': (k_rmse_kept + 1) / (N_ITERATIONS + 1)})
    for h in held_out:
        study_r2s, study_rmses = per[f'r2_{h}'], per[f'rmse_{h}']
        reach_r2, reach_rmse = study_r2s >= observed[h][0], study_rmses <= observed[h][1]
        k_r2_h, k_rmse_h = int(reach_r2.sum()), int(reach_rmse.sum())
        disagree[(scheme, h)] = int((reach_r2 != reach_rmse).sum())
        by_study.append({'shuffle': scheme, 'held_out_study': h, 'n_shuffles': N_ITERATIONS,
                         'observed_r2': observed[h][0], 'null_r2_mean': study_r2s.mean(),
                         'null_r2_median': study_r2s.median(),
                         'k_r2': k_r2_h, 'p_r2': (k_r2_h + 1) / (N_ITERATIONS + 1),
                         'observed_rmse': observed[h][1], 'null_rmse_mean': study_rmses.mean(),
                         'k_rmse': k_rmse_h, 'p_rmse': (k_rmse_h + 1) / (N_ITERATIONS + 1)})
    print(f"  Saved full results to: {NULL_PATHS[scheme]}")

print("-" * 80)
print(f"NULL DISTRIBUTIONS (N={N_ITERATIONS}), observed R2 {observed_r2:.3f}, RMSE {observed_rmse:.1f} mM C:")
for s in summary:
    print(f"  Shuffled {s['shuffle']:28s} R2 {s['null_r2_mean']:7.3f} (± {s['null_r2_sd']:.3f}), "
          f"median {s['null_r2_median']:7.3f} | k = {s['k_r2']:3d}, p = {s['p_r2']:.4f} || "
          f"RMSE {s['null_rmse_mean']:6.1f} (± {s['null_rmse_sd']:.1f}) | k = {s['k_rmse']:3d}, p = {s['p_rmse']:.4f}")
    print(f"    without the {LEFT_OUT_STUDY} held-out fold: k(R2) = {s['k_r2_without']}, p(R2) "
          f"{s['p_r2_without']:.4f} | k(RMSE) = {s['k_rmse_without']}, p(RMSE) {s['p_rmse_without']:.4f}")
print("  k = shuffles that do at least as well as the true labels (R2 as high or higher, or RMSE as low or lower)")
print("Per held-out study (observed R2 and RMSE; p of each shuffle scheme):")
for h in held_out:
    rows = [b for b in by_study if b['held_out_study'] == h]
    print(f"  {h:34s} R2 {observed[h][0]:7.3f}, RMSE {observed[h][1]:6.1f} | " +
          " | ".join(f"{b['shuffle']}: p(R2) {b['p_r2']:.4f}, p(RMSE) {b['p_rmse']:.4f}" for b in rows))
print("Self-check: per held-out study, p by R2 must equal p by RMSE when every sampling unit of the study "
      "has the same number of rows.")
if UNEVEN_STUDIES:
    print(f"  Studies whose sampling units differ in rows (no check): {', '.join(UNEVEN_STUDIES)}")
mismatches = [b for b in by_study
              if b['held_out_study'] not in UNEVEN_STUDIES and disagree[(b['shuffle'], b['held_out_study'])]]
for b in mismatches:
    print(f"  MISMATCH {b['shuffle']}, {b['held_out_study']}: {disagree[(b['shuffle'], b['held_out_study'])]} "
          f"shuffle(s) reach the observed R2 or RMSE but not both; k(R2) = {b['k_r2']}, k(RMSE) = {b['k_rmse']}")
n_checked = sum(1 for b in by_study if b['held_out_study'] not in UNEVEN_STUDIES)
print(f"  {n_checked - len(mismatches)} of {n_checked} study-and-shuffle pairs agree; {len(mismatches)} mismatch(es).")
pd.DataFrame(summary).to_csv(SUMMARY_PATH, index=False)
pd.DataFrame(by_study).to_csv(BY_STUDY_PATH, index=False)
print(f"  Saved: {SUMMARY_PATH}")
print(f"  Saved: {BY_STUDY_PATH}")
print("-" * 80)
