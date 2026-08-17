import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
import time
import os
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# CONFIG
# ==============================================================================
MATRIX_PATH = "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"
OUTPUT_DIR = "../../results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

df = pd.read_csv(MATRIX_PATH, low_memory=False)

genus_cols = sorted([col for col in df.columns if col.startswith('g__')])
df[genus_cols] = df[genus_cols].fillna(0)
row_sums = df[genus_cols].sum(axis=1).replace(0, 1)
df[genus_cols] = df[genus_cols].div(row_sums, axis=0)

TARGET = 'Caproate'
cat_cols = ['Feed_Complexity', 'Primary_Carbon_Signature']
for c in cat_cols:
    if c in df.columns:
        df[c] = df[c].astype('category')

papers_cstr = df['Paper_ID'].astype(str).str.strip()
y_true = pd.to_numeric(df[TARGET], errors='coerce').fillna(0)

all_chems = ['Lactate', 'Acetate', 'Ethanol', 'Butyrate', 'Valerate', 'Isovalerate', 'Propionate', 'Caprylate', 'Heptanoate']
base_ops = sorted([col for col in df.columns if col not in genus_cols + cat_cols + ['Paper_ID', 'Sample_ID', 'BIOREACTOR', 'Operation_Mode', 'DAY', TARGET, 'succinate'] + all_chems and not pd.api.types.is_string_dtype(df[col])])

valid_features = sorted(list(set(base_ops + cat_cols + genus_cols + all_chems)))
X = df[valid_features].copy()
numeric_cols = [col for col in X.columns if col not in cat_cols]
for c in numeric_cols:
    X[c] = pd.to_numeric(X[c], errors='coerce').fillna(0)

TRUE_OBSERVED_R2 = -0.095
N_ITERATIONS = 999  # <-- set to a small number (e.g. 20) to sanity-check first, then run the real 999

print(f"\n--- EXECUTING FORMAL {N_ITERATIONS}-ITERATION TARGET PERMUTATION TEST ---")
print(f"Benchmarking against Observed True-Label LOGO R2: {TRUE_OBSERVED_R2}")
print("-" * 80)

logo = LeaveOneGroupOut()
xgb_params = {'n_estimators': 150, 'learning_rate': 0.05, 'max_depth': 3, 'random_state': 42, 'n_jobs': -1, 'enable_categorical': True, 'tree_method': 'hist'}
permuted_r2_scores = []
start_time = time.time()

for iteration in range(1, N_ITERATIONS + 1):
    y_permuted = y_true.sample(frac=1.0, random_state=iteration).reset_index(drop=True)
    fold_r2 = []

    for train_idx, test_idx in logo.split(X, y_permuted, groups=papers_cstr):
        X_train, X_test = X.iloc[train_idx].copy(), X.iloc[test_idx].copy()
        y_train, y_test = y_permuted.iloc[train_idx], y_permuted.iloc[test_idx]

        scaler = StandardScaler()
        if len(numeric_cols) > 0:
            X_train[numeric_cols] = scaler.fit_transform(X_train[numeric_cols])
            X_test[numeric_cols] = scaler.transform(X_test[numeric_cols])

        model = xgb.XGBRegressor(**xgb_params).fit(X_train, y_train)
        y_pred = model.predict(X_test)
        if len(y_test) > 1:
            fold_r2.append(r2_score(y_test, y_pred))

    iter_mean_r2 = np.mean(fold_r2)
    permuted_r2_scores.append(iter_mean_r2)

    if iteration % 10 == 0 or iteration == 1 or iteration == N_ITERATIONS:
        elapsed = time.time() - start_time
        print(f"Iteration {iteration:03d}/{N_ITERATIONS} | Permuted LOGO R2: {iter_mean_r2:6.3f} | Elapsed Time: {elapsed:.1f}s")

# ==============================================================================
# NEW: save EVERY iteration's result (not just the printed subset) -- this is
# the file Figure 2 will load. This is also the file that belongs in the
# GitHub repository's results/ directory once we get there.
# ==============================================================================
results_df = pd.DataFrame({
    'iteration': range(1, N_ITERATIONS + 1),
    'permuted_r2': permuted_r2_scores,
})
save_path = os.path.join(OUTPUT_DIR, f'permutation_null_distribution_{N_ITERATIONS}.csv')
results_df.to_csv(save_path, index=False)

k_beats = sum(r2 >= TRUE_OBSERVED_R2 for r2 in permuted_r2_scores)
empirical_p_val = (k_beats + 1) / (N_ITERATIONS + 1)

print("-" * 80)
print(f"NULL DISTRIBUTION SUMMARY (N={N_ITERATIONS}):")
print(f"  Mean Permuted R2 : {np.mean(permuted_r2_scores):.3f} (\u00b1 {np.std(permuted_r2_scores):.3f})")
print(f"  True Observed R2 : {TRUE_OBSERVED_R2:.3f}")
print(f"  Shuffles Beating True Model (k) : {k_beats} / {N_ITERATIONS}")
print(f"  --> EMPIRICAL P-VALUE : p = {empirical_p_val:.4f}")
print(f"  Saved full results to: {save_path}")
print("-" * 80)
