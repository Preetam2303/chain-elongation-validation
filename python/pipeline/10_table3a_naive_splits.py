# 10_table3a_naive_splits.py
# Table 3a: the naive, ungrouped validation schemes, which the manuscript
# contrasts with the grouped schemes of Table 3b. Until 2026-10-10 these
# numbers came from code outside the repository; this script makes them.
#
#   1. Single random 80/20 split, ASV level (Illumina rows, BIOTWIN_PRUNED)
#   2. Random 5-fold CV, ASV level (same rows and predictors as 1)
#   3. Shuffled 10-fold CV with per-fold SHAP top-25 genera (full matrix,
#      BIOTWIN_MATRIX): the matched-feature-space comparator for Table 3b
#   4. Single-study internal 80/20 split (Duber 2024 rows of the full matrix)
#
# Every scheme uses the core predictor set (pH, temperature, HRT, lactate,
# acetate, ethanol, butyrate, propionate, the 2 categoricals) plus ASVs or
# genera, and the shared settings in pipeline_settings.py (no subsampling).
# Numeric predictors are scaled on the training rows only.
#
# Each scheme also reports how many same-day replicate sets (Duber 2024 and
# Nanopore triplicates) end up on both sides of a split -- the reason these
# scores are inflated (decision 2026-10-07).

import os
import warnings

import numpy as np
import pandas as pd
import shap
import xgboost as xgb
from sklearn.model_selection import KFold, train_test_split

from pipeline_settings import (CAT_COLS, CORE_CHEM, OPS, REPLICATE_STUDIES, TARGET, XGB_PARAMS, load_matrix,
                               present, r2_rmse, sampling_units, save_table, scale)

warnings.filterwarnings('ignore')

MATRIX_PATH = os.environ.get("BIOTWIN_MATRIX", "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")
PRUNED_PATH = os.environ.get("BIOTWIN_PRUNED", "../../data/historical/BIOTWIN_PRUNED_ML_MATRIX.csv")
CORE_NON_GENUS = OPS + CORE_CHEM + CAT_COLS


def load_asv_matrix(path):
    df = pd.read_csv(path, low_memory=False)
    for c in ['Paper_ID', 'Sample_ID', 'BIOREACTOR']:
        df[c] = df[c].astype(str).str.strip()
    for c in OPS + CORE_CHEM + [TARGET]:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    for c in CAT_COLS:
        df[c] = df[c].astype('category')
    if df[TARGET].isna().any():
        raise SystemExit("Caproate is missing in the ASV-level matrix; every modelled row needs it.")
    return df, [c for c in df.columns if c.startswith('ASV_')]


def split_sets(units, fold_of_row):
    """Replicate sets (units with more than one row) whose rows fall in more than one fold."""
    t = pd.DataFrame({'u': units.values, 'f': np.asarray(fold_of_row)})
    multi = t.groupby('u').filter(lambda g: len(g) > 1)
    if multi.empty:
        return 0, 0
    n_folds = multi.groupby('u')['f'].nunique()
    return int((n_folds > 1).sum()), int(len(n_folds))


def single_split(df, features, units):
    idx_train, idx_test = train_test_split(np.arange(len(df)), test_size=0.2, random_state=42)
    X_tr, X_te = scale(df[features].iloc[idx_train], df[features].iloc[idx_test])
    model = xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr, df[TARGET].iloc[idx_train])
    r2, rmse = r2_rmse(df[TARGET].iloc[idx_test], model.predict(X_te))
    side = np.zeros(len(df), int); side[idx_test] = 1
    split, total = split_sets(units, side)
    return {'r2_mean': r2, 'r2_sd': np.nan, 'rmse_mean': rmse, 'rmse_sd': np.nan, 'n': len(df),
            'n_test': len(idx_test), 'n_features': len(features),
            'replicate_sets_split': split, 'replicate_sets': total}


def kfold(df, features, units, k, top25=False):
    r2s, rmses = [], []
    fold_of_row = np.zeros(len(df), int)
    for f, (tr, te) in enumerate(KFold(k, shuffle=True, random_state=42).split(df), 1):
        fold_of_row[te] = f
        X_tr, X_te = scale(df[features].iloc[tr], df[features].iloc[te])
        cols = features
        if top25:
            explore = xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr, df[TARGET].iloc[tr])
            sv = np.abs(shap.TreeExplainer(explore).shap_values(X_tr)).mean(axis=0)
            ranked = pd.Series(sv, index=features).sort_values(ascending=False, kind='mergesort')
            top = [c for c in ranked.index if c.startswith('g__')][:25]
            cols = sorted(set([c for c in features if not c.startswith('g__')] + top))
        model = xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr[cols], df[TARGET].iloc[tr])
        r2, rmse = r2_rmse(df[TARGET].iloc[te], model.predict(X_te[cols]))
        r2s.append(r2); rmses.append(rmse)
    split, total = split_sets(units, fold_of_row)
    return {'r2_mean': np.mean(r2s), 'r2_sd': np.std(r2s), 'rmse_mean': np.mean(rmses), 'rmse_sd': np.std(rmses),
            'n': len(df), 'n_test': np.nan, 'n_features': len(features),
            'replicate_sets_split': split, 'replicate_sets': total,
            'folds': '; '.join(f"{r:.3f}" for r in r2s)}


asv, asv_cols = load_asv_matrix(PRUNED_PATH)
asv_units = sampling_units(asv)
asv_features = sorted(set(present(CORE_NON_GENUS, asv))) + asv_cols

full, genus_cols = load_matrix(MATRIX_PATH)
full_units = sampling_units(full)
full_features = sorted(set(present(CORE_NON_GENUS, full) + genus_cols))

d24 = full[full['Paper_ID'] == 'Duber_2024'].reset_index(drop=True)
d24_units = sampling_units(d24)

runs = [
    ("Single random 80/20 split (ASV level)", single_split(asv, asv_features, asv_units)),
    ("Random 5-fold CV (ASV level)", kfold(asv, asv_features, asv_units, 5)),
    ("Shuffled 10-fold CV, per-fold SHAP top-25 genera", kfold(full, full_features, full_units, 10, top25=True)),
    ("Single-study internal 80/20 split (Duber 2024)", single_split(d24, full_features, d24_units)),
]

print("--- TABLE 3a: NAIVE, UNGROUPED VALIDATION SCHEMES (core set, no subsampling) ---")
print(f"Replicate sets = same-day replicate samples of {', '.join(REPLICATE_STUDIES)}")
print("-" * 110)
table = []
for label, r in runs:
    sd = '' if np.isnan(r['r2_sd']) else f" (± {r['r2_sd']:.3f})"
    rsd = '' if np.isnan(r['rmse_sd']) else f" (± {r['rmse_sd']:.1f})"
    print(f"  {label:<52s} n {r['n']:3d} | {r['n_features']:5d} features | R2 {r['r2_mean']:6.3f}{sd} | "
          f"RMSE {r['rmse_mean']:6.1f}{rsd} mM C | replicate sets split: {r['replicate_sets_split']} of {r['replicate_sets']}")
    if 'folds' in r:
        print(f"      folds: {r['folds']}")
    table.append({'table': 'Table 3a', 'scheme': label, **r})
print("-" * 110)
save_table(table, "10_table3a_naive_splits.csv")
