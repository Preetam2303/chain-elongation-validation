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
# and all vessel types including batch reactors (Mode 2).

import pandas as pd
import numpy as np
import xgboost as xgb
import shap
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, root_mean_squared_error
import warnings
warnings.filterwarnings('ignore')

MATRIX_PATH = "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"
df = pd.read_csv(MATRIX_PATH)

genus_cols = [c for c in df.columns if c.startswith('g__')]
df[genus_cols] = df[genus_cols].fillna(0)
row_sums = df[genus_cols].sum(axis=1).replace(0, 1)
df[genus_cols] = df[genus_cols].div(row_sums, axis=0)

df['Unique_Vessel_ID'] = df['Paper_ID'].astype(str).str.strip() + "_" + df['BIOREACTOR'].astype(str).str.strip()

TARGET = 'Caproate'
foundational_chems = ['Lactate', 'Acetate', 'Butyrate', 'Ethanol', 'Propionate']
op_params = ['PH', 'TEMP', 'HRT']
valid_chems_ops = [c for c in foundational_chems + op_params if c in df.columns]

for c in valid_chems_ops + genus_cols + [TARGET]:
    df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0)

test_modes = {
    "Mode 1: CSTR Only (B1, B2 Vessels)": df[df['BIOREACTOR'].astype(str).str.strip().str.startswith('B')].copy(),
    "Mode 2: All Vessels (CSTR + Batch R1-R9)": df.copy(),
}

print("--- EXECUTING VESSEL-GROUPED DEPLOYABILITY TEST (Table 3b) ---")
print("-" * 85)

xgb_params = {
    'n_estimators': 150, 'learning_rate': 0.05, 'max_depth': 3,
    'subsample': 0.8, 'colsample_bytree': 0.8, 'random_state': 42,
    'n_jobs': -1, 'tree_method': 'hist',
}

for mode_name, subset_df in test_modes.items():
    print(f"\n{mode_name}")
    print(f"Total Samples: {len(subset_df)} | Unique Physical Vessels: {subset_df['Unique_Vessel_ID'].nunique()}")

    X = subset_df[valid_chems_ops + genus_cols].copy()
    y = subset_df[TARGET].copy()
    groups = subset_df['Unique_Vessel_ID'].copy()

    n_splits = min(5, subset_df['Unique_Vessel_ID'].nunique())
    gkf = GroupKFold(n_splits=n_splits)

    r2_scores, rmse_scores = [], []
    fold = 1

    for train_idx, test_idx in gkf.split(X, y, groups=groups):
        X_train_full, X_test_full = X.iloc[train_idx].copy(), X.iloc[test_idx].copy()
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        test_vessels = groups.iloc[test_idx].unique().tolist()

        # Train-only scaling
        scaler = StandardScaler()
        X_train_full[valid_chems_ops + genus_cols] = scaler.fit_transform(X_train_full[valid_chems_ops + genus_cols])
        X_test_full[valid_chems_ops + genus_cols] = scaler.transform(X_test_full[valid_chems_ops + genus_cols])

        # Phase A: dynamic SHAP feature discovery on training vessels only
        model_discovery = xgb.XGBRegressor(**xgb_params).fit(X_train_full, y_train)
        explainer = shap.TreeExplainer(model_discovery)
        shap_values = explainer.shap_values(X_train_full)

        mean_shap = np.abs(shap_values).mean(axis=0)
        shap_df = pd.DataFrame({'Feature': X_train_full.columns, 'SHAP': mean_shap})
        top_25_genera = shap_df.sort_values(by='SHAP', ascending=False)[shap_df['Feature'].str.startswith('g__')].head(25)['Feature'].tolist()

        # Phase B: final prediction on held-out physical vessels
        final_features = valid_chems_ops + top_25_genera
        model_prediction = xgb.XGBRegressor(**xgb_params).fit(X_train_full[final_features], y_train)
        y_pred = model_prediction.predict(X_test_full[final_features])

        if len(y_test) > 1:
            r2 = r2_score(y_test, y_pred)
            rmse = root_mean_squared_error(y_test, y_pred)
            r2_scores.append(r2)
            rmse_scores.append(rmse)
            vessel_str = ", ".join(test_vessels[:3]) + ("..." if len(test_vessels) > 3 else "")
            print(f"  Fold {fold} | Held-Out Vessels: {vessel_str:<30} | R2: {r2:6.3f} | RMSE: {rmse:5.1f} mM C")
        fold += 1

    print(f"  --> RESULT: Average R2: {np.mean(r2_scores):.3f} (\u00b1 {np.std(r2_scores):.3f}) | "
          f"Average RMSE: {np.mean(rmse_scores):.1f} mM C")
    print("-" * 85)
