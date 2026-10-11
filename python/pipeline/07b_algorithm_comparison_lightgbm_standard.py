# 07b_algorithm_comparison_lightgbm_standard.py
# The second half of the algorithm-selection comparison in Methods 2.7: a
# standard gradient-boosted (gbdt) LightGBM configuration against XGBoost.
# Both models use the shared settings (150 trees, learning rate 0.05, depth 3,
# no subsampling), so the comparison is like for like.
#
# The same two models in three validation schemes (decision 2026-10-11: the
# ranking of the two algorithms depends on the scheme):
#   1. Single 80/20 split, ASV level: the exact feature set and split
#      (random_state=42) of the Table 3a single-split baseline
#      (10_table3a_naive_splits.py), on BIOTWIN_PRUNED.
#   2. Shuffled 10-fold CV and
#   3. leave-one-study-out (LOSO), both with the headline's rows and predictors
#      (02_loso_by_study.py: full matrix, core set + genera as relative
#      abundance, predictors scaled on the training rows only), on
#      BIOTWIN_MATRIX. Only the way rows are split differs between 2 and 3;
#      XGBoost's LOSO is the headline itself. The training-mean benchmark is
#      shown with both.
#   4. The two naive schemes (1 and 2) again with seeds 0-19. Which rows land
#      in the test set changes a naive score a lot, so one seed can rank the
#      two algorithms either way; LOSO has no random split.
#
# Feature set: the core set (pH, temperature, HRT, lactate, acetate, ethanol,
# butyrate, propionate, the 2 categoricals) + ASVs or genera, as everywhere
# else (pipeline_settings.py). Unmeasured cells stay blank; both libraries
# handle blanks.
#
# Output (in BIOTWIN_RESULTS_DIR): 07b_algorithm_comparison.csv, one row per
# scheme and model, with each fold's R2, and
# 07b_algorithm_comparison_over_seeds.csv, one row per seed, scheme and model.
#
# Published before this change: XGBoost 0.915 (typed in, not computed),
# LightGBM 0.861, from a feature set that included NaOH and DAY.

import os
import warnings

import numpy as np
import pandas as pd
import xgboost as xgb
from lightgbm import LGBMRegressor
from sklearn.metrics import r2_score, root_mean_squared_error
from sklearn.model_selection import KFold, LeaveOneGroupOut, train_test_split

from pipeline_settings import (CAT_COLS, CORE_CHEM, OPS, TARGET, XGB_PARAMS, drop_unmodelled, load_matrix,
                               present, r2_rmse, save_table, scale, training_mean_r2_rmse)

warnings.filterwarnings('ignore')

# 'gbdt' (Gradient Boosting Decision Tree), not 'rf' -- this is the standard
# boosting comparison, distinct from 07a's bagged Random-Forest-mode run.
LGBM_PARAMS = dict(
    boosting_type='gbdt', n_estimators=150, learning_rate=0.05, max_depth=3,
    num_leaves=7,             # capped safely below 2^max_depth
    min_child_samples=5,      # minimum samples in a leaf to allow a split
    random_state=42, n_jobs=-1, verbose=-1,
)

# ==============================================================================
# 1. LOAD THE PRUNED MATRIX & ENVIRONMENT SETUP
# ==============================================================================
file_path = os.environ.get("BIOTWIN_PRUNED", "../../data/historical/BIOTWIN_PRUNED_ML_MATRIX.csv")
print("Loading Pruned ML Matrix...")
df = drop_unmodelled(pd.read_csv(file_path, low_memory=False))

core_asvs = [col for col in df.columns if str(col).startswith('ASV_')]

# ==============================================================================
# 2. FEATURE SET & TARGET
# ==============================================================================
core_cols = OPS + CORE_CHEM + CAT_COLS
for col in OPS + CORE_CHEM + [TARGET]:
    df[col] = pd.to_numeric(df[col], errors='coerce')
for col in CAT_COLS:
    df[col] = df[col].astype('category')
if df[TARGET].isna().any():
    raise SystemExit("Caproate is missing in the ASV-level matrix; every modelled row needs it.")

# The same columns, in the same order, as the Table 3a single-split baseline
X = df[sorted(set(present(core_cols, df))) + core_asvs]
y = df[TARGET]
# The same 80/20 split as the Table 3a single-split baseline
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

print(f"Training Matrix built with {X_train.shape[1]} features.")
print(f"Training on {X_train.shape[0]} samples, testing on {X_test.shape[0]} samples.")

# ==============================================================================
# 3. XGBOOST (shared settings) AND STANDARD GRADIENT-BOOSTED LIGHTGBM
# ==============================================================================
xgb_pred = xgb.XGBRegressor(**XGB_PARAMS).fit(X_train, y_train).predict(X_test)

print("\nTraining Standard LightGBM (GBDT) on the core set + ASVs...")
lgbm_pred = LGBMRegressor(**LGBM_PARAMS).fit(X_train, y_train).predict(X_test)

# ==============================================================================
# 4. RESULTS
# ==============================================================================
print("\n" + "=" * 50)
print("XGBOOST VS STANDARD LIGHTGBM, SAME FEATURES AND SPLIT")
print("=" * 50)
print(f"Target Variable:     {TARGET}")
print(f"XGBoost R-squared:   {r2_score(y_test, xgb_pred):.3f}")
print(f"XGBoost RMSE:        {root_mean_squared_error(y_test, xgb_pred):.3f} mM C")
print(f"LightGBM R-squared:  {r2_score(y_test, lgbm_pred):.3f}")
print(f"LightGBM RMSE:       {root_mean_squared_error(y_test, lgbm_pred):.3f} mM C")
print("=" * 50)

split_label = "Single 80/20 split (ASV level)"
table = [{'scheme': split_label, 'model': model, 'n_rows': len(df), 'n_features': X.shape[1],
          'r2_mean': r2_score(y_test, pred), 'r2_sd': np.nan,
          'rmse_mean': root_mean_squared_error(y_test, pred), 'rmse_sd': np.nan, 'folds': ''}
         for model, pred in [('XGBoost', xgb_pred), ('LightGBM', lgbm_pred)]]

# ==============================================================================
# 5. THE HEADLINE'S ROWS AND PREDICTORS: SHUFFLED 10-FOLD VS LOSO
# ==============================================================================
MATRIX_PATH = os.environ.get("BIOTWIN_MATRIX", "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")
full, genus_cols = load_matrix(MATRIX_PATH)
features = sorted(set(present(OPS + CORE_CHEM + CAT_COLS, full) + genus_cols))
X_full, y_full, papers = full[features], full[TARGET], full['Paper_ID']
MODELS = ['XGBoost', 'LightGBM', 'Training-mean benchmark']


def cross_validate(splits):
    """R2 and RMSE of each model in each fold. splits: (fold label, train rows, test rows)."""
    res = {m: {'r2': [], 'rmse': []} for m in MODELS}
    labels = []
    for label, tr, te in splits:
        X_tr, X_te = scale(X_full.iloc[tr], X_full.iloc[te])
        y_tr, y_te = y_full.iloc[tr], y_full.iloc[te]
        scores = {
            'XGBoost': r2_rmse(y_te, xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr, y_tr).predict(X_te)),
            'LightGBM': r2_rmse(y_te, LGBMRegressor(**LGBM_PARAMS).fit(X_tr, y_tr).predict(X_te)),
            'Training-mean benchmark': training_mean_r2_rmse(y_tr, y_te),
        }
        for m, (r2, rmse) in scores.items():
            res[m]['r2'].append(r2); res[m]['rmse'].append(rmse)
        labels.append(label)
    return labels, res


schemes = [
    ("Shuffled 10-fold CV (headline predictors)",
     [(f"fold {i}", tr, te) for i, (tr, te) in enumerate(KFold(10, shuffle=True, random_state=42).split(X_full), 1)]),
    ("LOSO (headline)",
     [(papers.iloc[te[0]], tr, te) for tr, te in LeaveOneGroupOut().split(X_full, y_full, groups=papers)
      if len(te) >= 2]),
]

print(f"\n--- SAME ROWS AND PREDICTORS AS THE HEADLINE: {len(full)} rows, {papers.nunique()} studies, "
      f"{len(features)} predictors ---")
for scheme, splits in schemes:
    labels, res = cross_validate(splits)
    print(f"\n{scheme}")
    for m in MODELS:
        r2s, rmses = res[m]['r2'], res[m]['rmse']
        print(f"  {m:24s} R2 {np.mean(r2s):7.3f} (± {np.std(r2s):.3f})  "
              f"RMSE {np.mean(rmses):6.1f} (± {np.std(rmses):.1f}) mM C")
        table.append({'scheme': scheme, 'model': m, 'n_rows': len(full),
                      'n_features': 0 if m == 'Training-mean benchmark' else len(features),
                      'r2_mean': np.mean(r2s), 'r2_sd': np.std(r2s),
                      'rmse_mean': np.mean(rmses), 'rmse_sd': np.std(rmses),
                      'folds': '; '.join(f"{lab} {r:.3f}" if scheme.startswith('LOSO') else f"{r:.3f}"
                                         for lab, r in zip(labels, r2s))})
    if scheme.startswith('LOSO'):
        print("  R2 per held-out study:")
        for i, lab in enumerate(labels):
            print(f"    {lab:34s} " + "  ".join(f"{m.split()[0]} {res[m]['r2'][i]:7.3f}" for m in MODELS))

print("\nXGBoost's LOSO row is the headline (Table 3b = Table 4 full-matrix tier 3, 02_loso_by_study.py).")
save_table(table, "07b_algorithm_comparison.csv")

# ==============================================================================
# 6. THE TWO NAIVE SCHEMES WITH SEEDS 0-19
# ==============================================================================
SEEDS = range(20)
over_seeds = []
for seed in SEEDS:
    X_tr_s, X_te_s, y_tr_s, y_te_s = train_test_split(X, y, test_size=0.2, random_state=seed)
    for model, pred in [('XGBoost', xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr_s, y_tr_s).predict(X_te_s)),
                        ('LightGBM', LGBMRegressor(**LGBM_PARAMS).fit(X_tr_s, y_tr_s).predict(X_te_s))]:
        over_seeds.append({'seed': seed, 'scheme': split_label, 'model': model, 'r2': r2_score(y_te_s, pred)})
    _, res = cross_validate([(f"fold {i}", tr, te) for i, (tr, te) in
                             enumerate(KFold(10, shuffle=True, random_state=seed).split(X_full), 1)])
    for model in ['XGBoost', 'LightGBM']:
        over_seeds.append({'seed': seed, 'scheme': schemes[0][0], 'model': model, 'r2': np.mean(res[model]['r2'])})

print(f"\n--- THE TWO NAIVE SCHEMES WITH SEEDS {min(SEEDS)}-{max(SEEDS)} (above: seed 42): "
      "R2 median [lowest, highest] ---")
seeds_df = pd.DataFrame(over_seeds)
for scheme, g in seeds_df.groupby('scheme', sort=False):
    w = g.pivot(index='seed', columns='model', values='r2')
    print(f"  {scheme}")
    for model in ['XGBoost', 'LightGBM']:
        print(f"    {model:9s} {w[model].median():6.3f} [{w[model].min():.3f}, {w[model].max():.3f}]")
    print(f"    LightGBM higher in {int((w['LightGBM'] > w['XGBoost']).sum())} of {len(w)} seeds; "
          f"median difference (LightGBM - XGBoost) {(w['LightGBM'] - w['XGBoost']).median():+.3f}")
save_table(over_seeds, "07b_algorithm_comparison_over_seeds.csv")
