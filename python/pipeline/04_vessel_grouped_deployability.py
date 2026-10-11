# 04_vessel_grouped_deployability.py
# GroupKFold by physical bioreactor vessel (Methods 2.6, scheme iii). Tests a
# narrower, practically motivated question than 02_loso_by_study.py: not "can
# this generalize to a lab never seen before" but "can periodic microbiome
# sampling substitute for continuous chemical monitoring within an
# already-characterized dataset." Groups by Study x Bioreactor ID so that all
# longitudinal timepoints from one physical tank -- including replicate rows
# from the same reactor-day -- are quarantined to a single fold, precluding
# both cross-study and within-reactor temporal leakage.
#
# Produces both rows of Table 3b's vessel-grouped section: CSTR-only (Mode 1)
# and all vessel types including batch reactors (Mode 2), each with the
# training-mean benchmark next to it (review item 5).
#
# In each fold, an exploratory model on the training vessels picks the 25
# genera with the largest mean |SHAP|; the final model uses the core chemistry,
# operational and categorical predictors plus those 25 genera. Settings as in
# pipeline_settings.py (no subsampling), so every number reproduces exactly.

import os
import warnings

import numpy as np
import pandas as pd
import shap
import xgboost as xgb
from sklearn.model_selection import GroupKFold

from pipeline_settings import (CAT_COLS, CORE_CHEM, OPS, TARGET, XGB_PARAMS, load_matrix, present,
                               r2_rmse, save_table, scale, training_mean_r2_rmse)

warnings.filterwarnings('ignore')

MATRIX_PATH = os.environ.get("BIOTWIN_MATRIX", "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")

df, genus_cols = load_matrix(MATRIX_PATH)
df['Unique_Vessel_ID'] = df['Paper_ID'] + "_" + df['BIOREACTOR']
core_non_genus = present(OPS + CORE_CHEM + CAT_COLS, df)
features = sorted(set(core_non_genus + genus_cols))

test_modes = {
    "Mode 1: CSTR Only (B1, B2 Vessels)": df[df['BIOREACTOR'].str.startswith('B')].reset_index(drop=True),
    "Mode 2: All Vessels (CSTR + Batch R1-R9)": df,
}

print("--- EXECUTING VESSEL-GROUPED DEPLOYABILITY TEST (Table 3b) ---")
print("-" * 85)

table = []
for mode_name, subset_df in test_modes.items():
    print(f"\n{mode_name}")
    print(f"Total Samples: {len(subset_df)} | Unique Physical Vessels: {subset_df['Unique_Vessel_ID'].nunique()}")

    X = subset_df[features]
    y = subset_df[TARGET]
    groups = subset_df['Unique_Vessel_ID']
    gkf = GroupKFold(n_splits=min(5, groups.nunique()))

    r2_scores, rmse_scores, base_r2, base_rmse = [], [], [], []
    for fold, (train_idx, test_idx) in enumerate(gkf.split(X, y, groups=groups), 1):
        X_train, X_test = scale(X.iloc[train_idx], X.iloc[test_idx])
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        test_vessels = groups.iloc[test_idx].unique().tolist()

        # Phase A: dynamic SHAP feature discovery on training vessels only
        discovery = xgb.XGBRegressor(**XGB_PARAMS).fit(X_train, y_train)
        mean_shap = np.abs(shap.TreeExplainer(discovery).shap_values(X_train)).mean(axis=0)
        ranked = pd.Series(mean_shap, index=features).sort_values(ascending=False, kind='mergesort')
        top_25_genera = [c for c in ranked.index if c.startswith('g__')][:25]

        # Phase B: final prediction on held-out physical vessels
        final_features = sorted(set(core_non_genus + top_25_genera))
        model = xgb.XGBRegressor(**XGB_PARAMS).fit(X_train[final_features], y_train)
        y_pred = model.predict(X_test[final_features])

        if len(y_test) > 1:
            r2, rmse = r2_rmse(y_test, y_pred)
            b_r2, b_rmse = training_mean_r2_rmse(y_train, y_test)
            r2_scores.append(r2); rmse_scores.append(rmse)
            base_r2.append(b_r2); base_rmse.append(b_rmse)
            vessel_str = ", ".join(test_vessels[:3]) + ("..." if len(test_vessels) > 3 else "")
            print(f"  Fold {fold} | Held-Out Vessels: {vessel_str:<30} | R2: {r2:6.3f} | RMSE: {rmse:5.1f} mM C"
                  f" | training-mean R2: {b_r2:7.3f}")

    print(f"  --> RESULT: Average R2: {np.mean(r2_scores):.3f} (± {np.std(r2_scores):.3f}) | "
          f"Average RMSE: {np.mean(rmse_scores):.1f} (± {np.std(rmse_scores):.1f}) mM C")
    print(f"  --> TRAINING-MEAN BENCHMARK: Average R2: {np.mean(base_r2):.3f} (± {np.std(base_r2):.3f}) | "
          f"Average RMSE: {np.mean(base_rmse):.1f} (± {np.std(base_rmse):.1f}) mM C")
    print("-" * 85)
    for model_name, r2s, rmses in [("XGBoost, core set + per-fold SHAP top-25 genera", r2_scores, rmse_scores),
                                   ("Training-mean benchmark", base_r2, base_rmse)]:
        table.append({'table': 'Table 3b', 'scheme': mode_name, 'model': model_name,
                      'n_rows': len(subset_df), 'n_vessels': groups.nunique(),
                      'r2_mean': np.mean(r2s), 'r2_sd': np.std(r2s),
                      'rmse_mean': np.mean(rmses), 'rmse_sd': np.std(rmses)})

save_table(table, "04_groupkfold_by_vessel.csv")
