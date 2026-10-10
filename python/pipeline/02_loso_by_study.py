# 02_loso_by_study.py
# Leave-One-Study-Out cross-validation (Methods 2.6, scheme ii): Table 4 (both
# columns), Table 3b's two LOSO rows, the chemistry-only LOSO, and the
# training-mean benchmark for each (review item 5).
#
# Full-matrix column (all 7 groups, Illumina + Nanopore). Only compounds
# measured in every study, so it stops at tier 3:
#   Tier 1  pH, temperature, HRT, the 2 categoricals, genera
#   Tier 2  + lactate, acetate, ethanol
#   Tier 3  + butyrate, propionate = exactly the core set = THE HEADLINE
#           (Table 3b's LOSO row shows this same number)
# Illumina-only column (the 6 Illumina studies):
#   Tiers 1-2 as above
#   Tier 3  + butyrate, propionate, valerate, isovalerate
#   Tier 4  + caprylate, heptanoate, isocaproate (the "even with downstream
#           proxies" check)
# Also:
#   - chemistry only: the core set without genera, full matrix;
#   - Table 3b "LOSO, SHAP top-25": in each fold, an exploratory model on the
#     training studies picks the 25 genera with the largest mean |SHAP|, and the
#     final model uses the core chemistry and operational predictors plus those
#     25 genera;
#   - training-mean benchmark: predict the training studies' mean caproate.
#     R2 is computed within each held-out study, against that study's own mean.
#
# Settings (pipeline_settings.py): relative abundance, no subsampling, numeric
# predictors scaled on the training studies only, target not scaled. Every
# number here reproduces exactly on any machine.
#
# Feature lists are built via sorted(set(...)) deliberately: Python's default
# set iteration order is process-randomized, and column order can change which
# of two equally good splits XGBoost picks.

import os
import warnings

import numpy as np
import pandas as pd
import shap
import xgboost as xgb
from sklearn.model_selection import LeaveOneGroupOut

from pipeline_settings import (CORE_CHEM, CAT_COLS, FEEDS, ILLUMINA_ONLY_TIER3, ILLUMINA_ONLY_TIER4,
                               INTERMEDIATES, NANOPORE_ID, OPS, TARGET, XGB_PARAMS, load_matrix,
                               present, r2_rmse, save_table, scale, training_mean_r2_rmse)

warnings.filterwarnings('ignore')

MATRIX_PATH = os.environ.get("BIOTWIN_MATRIX", "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")

df, genus_cols = load_matrix(MATRIX_PATH)
illumina = df[df['Paper_ID'] != NANOPORE_ID].reset_index(drop=True)


def loso(data, features, top25=False):
    """Mean and SD over held-out studies of R2 and RMSE, for the model and for
    the training-mean benchmark."""
    features = sorted(set(present(features, data)))
    X, y, groups = data[features], data[TARGET], data['Paper_ID']
    out = {'model_r2': [], 'model_rmse': [], 'mean_r2': [], 'mean_rmse': [], 'held_out': []}
    for tr, te in LeaveOneGroupOut().split(X, y, groups=groups):
        if len(te) < 2:
            continue
        X_tr, X_te = scale(X.iloc[tr], X.iloc[te])
        cols = features
        if top25:
            explore = xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr, y.iloc[tr])
            sv = np.abs(shap.TreeExplainer(explore).shap_values(X_tr)).mean(axis=0)
            ranked = pd.Series(sv, index=features).sort_values(ascending=False, kind='mergesort')
            top = [c for c in ranked.index if c.startswith('g__')][:25]
            cols = sorted(set([c for c in features if not c.startswith('g__')] + top))
        model = xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr[cols], y.iloc[tr])
        r2, rmse = r2_rmse(y.iloc[te], model.predict(X_te[cols]))
        b_r2, b_rmse = training_mean_r2_rmse(y.iloc[tr], y.iloc[te])
        out['model_r2'].append(r2); out['model_rmse'].append(rmse)
        out['mean_r2'].append(b_r2); out['mean_rmse'].append(b_rmse)
        out['held_out'].append(groups.iloc[te].iloc[0])
    return features, out


def line(label, vals_r2, vals_rmse):
    return (f"  {label:62s} R2 {np.mean(vals_r2):7.3f} (± {np.std(vals_r2):.3f})  "
            f"RMSE {np.mean(vals_rmse):6.1f} (± {np.std(vals_rmse):.1f}) mM C")


tier1 = OPS + CAT_COLS + genus_cols
tier2 = tier1 + FEEDS
tier3_core = tier2 + INTERMEDIATES
runs = [
    # (Table, column, label, data, features, top25)
    ("Table 4", "full matrix", "Tier 1: pH, temperature, HRT, categoricals + genera", df, tier1, False),
    ("Table 4", "full matrix", "Tier 2: + lactate, acetate, ethanol", df, tier2, False),
    ("Table 4", "full matrix", "Tier 3: + butyrate, propionate (core set; HEADLINE)", df, tier3_core, False),
    ("Table 4", "Illumina-only", "Tier 1: pH, temperature, HRT, categoricals + genera", illumina, tier1, False),
    ("Table 4", "Illumina-only", "Tier 2: + lactate, acetate, ethanol", illumina, tier2, False),
    ("Table 4", "Illumina-only", "Tier 3: + butyrate, propionate, valerate, isovalerate", illumina,
     tier3_core + ILLUMINA_ONLY_TIER3, False),
    ("Table 4", "Illumina-only", "Tier 4: + caprylate, heptanoate, isocaproate", illumina,
     tier3_core + ILLUMINA_ONLY_TIER3 + ILLUMINA_ONLY_TIER4, False),
    ("Table 4", "full matrix", "Chemistry only: core set without genera", df, OPS + CORE_CHEM + CAT_COLS, False),
    ("Table 3b", "full matrix", "LOSO, core set, per-fold SHAP top-25 genera", df, tier3_core, True),
]

print("--- LEAVE-ONE-STUDY-OUT (Table 4, Table 3b LOSO rows) ---")
print(f"Full matrix: {len(df)} rows, {df['Paper_ID'].nunique()} studies, {len(genus_cols)} genera")
print(f"Illumina-only: {len(illumina)} rows, {illumina['Paper_ID'].nunique()} studies")
print("-" * 110)

table = []
benchmark_done = set()
previous = None
for table_name, column, label, data, features, top25 in runs:
    used, res = loso(data, features, top25)
    if column != previous:
        print(f"\n{column} ({len(data)} rows, {data['Paper_ID'].nunique()} held-out studies)")
        previous = column
    if column not in benchmark_done:
        print(line("Training-mean benchmark (predict the training studies' mean)", res['mean_r2'], res['mean_rmse']))
        table.append({'table': 'Tables 3b and 4', 'column': column, 'model': 'Training-mean benchmark',
                      'n_rows': len(data), 'n_features': 0,
                      'r2_mean': np.mean(res['mean_r2']), 'r2_sd': np.std(res['mean_r2']),
                      'rmse_mean': np.mean(res['mean_rmse']), 'rmse_sd': np.std(res['mean_rmse']),
                      'folds': '; '.join(f"{h} {r:.3f}" for h, r in zip(res['held_out'], res['mean_r2']))})
        benchmark_done.add(column)
    print(line(label, res['model_r2'], res['model_rmse']))
    table.append({'table': table_name, 'column': column, 'model': label,
                  'n_rows': len(data), 'n_features': len(used),
                  'r2_mean': np.mean(res['model_r2']), 'r2_sd': np.std(res['model_r2']),
                  'rmse_mean': np.mean(res['model_rmse']), 'rmse_sd': np.std(res['model_rmse']),
                  'folds': '; '.join(f"{h} {r:.3f}" for h, r in zip(res['held_out'], res['model_r2']))})

print("-" * 110)
headline = next(t for t in table if 'HEADLINE' in t['model'])
print(f"HEADLINE (full-matrix LOSO, core set; Table 3b = Table 4 full-matrix tier 3): "
      f"R2 = {headline['r2_mean']:.3f} ± {headline['r2_sd']:.3f}, "
      f"RMSE = {headline['rmse_mean']:.1f} ± {headline['rmse_sd']:.1f} mM C")
print(f"  per held-out study: {headline['folds']}")
save_table(table, "02_loso_tables_3b_4.csv")
