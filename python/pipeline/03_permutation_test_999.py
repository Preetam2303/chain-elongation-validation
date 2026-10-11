# 03_permutation_test_999.py
# Target-permutation control for the headline (Section 3.2, Figure 2): the
# full-matrix LOSO with exactly the core predictor set, the same model as
# Table 3b's LOSO row and Table 4's full-matrix tier 3 (02_loso_by_study.py).
#
# Caproate is shuffled by sampling unit, not by row (decision 2026-10-08): the
# three same-day replicate samples of Duber 2024 and of the Nanopore study
# share one caproate value, so they move together; every other row is its own
# unit. Shuffle i uses random_state = i, as before. The observed value is
# computed in the same run, with the true labels and the same features.
#
# p = (k + 1) / (N + 1), where k is the number of shuffles whose mean LOSO R2
# reaches the observed one.
#
# Output (in BIOTWIN_RESULTS_DIR, default ../../results):
#   permutation_null_distribution_<N>.csv   every shuffle's mean LOSO R2
# Set BIOTWIN_N_PERMUTATIONS (default 999) to e.g. 20 for a quick check.

import os
import time
import warnings

import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.metrics import r2_score
from sklearn.model_selection import LeaveOneGroupOut

from pipeline_settings import (CAT_COLS, CORE_CHEM, OPS, TARGET, XGB_PARAMS, load_matrix, present,
                               sampling_units, scale)

warnings.filterwarnings('ignore')

MATRIX_PATH = os.environ.get("BIOTWIN_MATRIX", "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")
OUTPUT_DIR = os.environ.get("BIOTWIN_RESULTS_DIR", "../../results")
os.makedirs(OUTPUT_DIR, exist_ok=True)
N_ITERATIONS = int(os.environ.get("BIOTWIN_N_PERMUTATIONS", 999))
save_path = os.path.join(OUTPUT_DIR, f'permutation_null_distribution_{N_ITERATIONS}.csv')
if os.path.exists(save_path):  # checked before the long run, not after it
    raise SystemExit(f"{save_path} already exists. Set BIOTWIN_RESULTS_DIR to a new folder.")

df, genus_cols = load_matrix(MATRIX_PATH)
features = sorted(set(present(OPS + CORE_CHEM + CAT_COLS, df) + genus_cols))
X = df[features]
papers = df['Paper_ID']
y_true = df[TARGET]

# Precompute the folds and their scaled predictors once; only y changes per shuffle.
folds = []
for train_idx, test_idx in LeaveOneGroupOut().split(X, y_true, groups=papers):
    if len(test_idx) > 1:
        X_train, X_test = scale(X.iloc[train_idx], X.iloc[test_idx])
        folds.append((train_idx, test_idx, X_train, X_test))


def loso_mean_r2(y):
    fold_r2 = []
    for train_idx, test_idx, X_train, X_test in folds:
        model = xgb.XGBRegressor(**XGB_PARAMS).fit(X_train, y.iloc[train_idx])
        fold_r2.append(r2_score(y.iloc[test_idx], model.predict(X_test)))
    return np.mean(fold_r2)


units = sampling_units(df)
unit_values = y_true.groupby(units).first()   # one caproate value per sampling unit
unit_order = unit_values.index


def shuffled_target(seed):
    shuffled = pd.Series(unit_values.sample(frac=1.0, random_state=seed).values, index=unit_order)
    return pd.Series(units.map(shuffled).values, index=df.index)


TRUE_OBSERVED_R2 = round(loso_mean_r2(y_true), 3)

print(f"\n--- {N_ITERATIONS}-SHUFFLE TARGET PERMUTATION TEST (headline: full-matrix LOSO, core set) ---")
print(f"Rows: {len(df)} | sampling units shuffled: {len(unit_values)} "
      f"(same-day replicates move together) | predictors: {len(features)}")
print(f"Observed LOSO R2 (true labels): {TRUE_OBSERVED_R2}")
print("-" * 80)

permuted_r2_scores = []
start_time = time.time()
for iteration in range(1, N_ITERATIONS + 1):
    iter_mean_r2 = loso_mean_r2(shuffled_target(iteration))
    permuted_r2_scores.append(iter_mean_r2)
    if iteration % 10 == 0 or iteration == 1 or iteration == N_ITERATIONS:
        elapsed = time.time() - start_time
        print(f"Iteration {iteration:03d}/{N_ITERATIONS} | Permuted LOGO R2: {iter_mean_r2:6.3f} | Elapsed Time: {elapsed:.1f}s")

results_df = pd.DataFrame({'iteration': range(1, N_ITERATIONS + 1), 'permuted_r2': permuted_r2_scores})
results_df.to_csv(save_path, index=False)

k_beats = sum(r2 >= TRUE_OBSERVED_R2 for r2 in permuted_r2_scores)
empirical_p_val = (k_beats + 1) / (N_ITERATIONS + 1)

print("-" * 80)
print(f"NULL DISTRIBUTION SUMMARY (N={N_ITERATIONS}):")
print(f"  Mean Permuted R2 : {np.mean(permuted_r2_scores):.3f} (± {np.std(permuted_r2_scores):.3f}); "
      f"median {np.median(permuted_r2_scores):.3f}")
print(f"  True Observed R2 : {TRUE_OBSERVED_R2:.3f}")
print(f"  Shuffles Beating True Model (k) : {k_beats} / {N_ITERATIONS}")
print(f"  --> EMPIRICAL P-VALUE : p = {empirical_p_val:.4f}")
print(f"  Saved full results to: {save_path}")
print("-" * 80)
