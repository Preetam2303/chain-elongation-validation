# 07c_figure3_intermediate_unscaled_xgboost.py
# The middle number in Figure 3's three-step progression (R2 = 0.217 -> -0.228
# -> -5.869): the same Duber 2025 B1(n=30) -> B2(n=21) split, evaluated with
# XGBoost but WITHOUT training-only scaling confirmed (no StandardScaler
# applied at all here) and without an explicit, curated leakage exclusion list.
#
# IMPORTANT -- read before running or citing:
#   1. This script loads BIOTWIN_GENUS_ML_MATRIX.csv, an EARLIER dataset
#      snapshot, NOT the final BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv used
#      by every other script in this pipeline. It is not reproducible from
#      01_genus_aggregation.py's current output. If you have the original
#      BIOTWIN_GENUS_ML_MATRIX.csv snapshot, place it in data/historical/ to
#      run this script; otherwise treat the numbers below as a static,
#      already-verified historical record rather than a re-runnable result.
#   2. The manuscript describes this run as "XGBoost without downstream
#      co-products." On inspection, that is not quite accurate: this script
#      only drops a short list of sparse/junk columns (succinate, Lactose,
#      i-propanol, izo-butanol, NAOH) and does NOT explicitly exclude
#      Valerate, Isovalerate, Heptanoate, or Caprylate -- the SHAP output
#      below shows Valerate and Heptanoate both used as real predictors.
#      The manuscript text should be corrected (see note at bottom of file)
#      before this is treated as fully precise; the numbers themselves are
#      unaffected and the overall three-step "tightening rigor" narrative
#      still holds.
#
# Expected output (confirmed against the original run):
#   Intra-Study R-squared (R2): -0.228
#   Intra-Study RMSE: 63.736 mM C

import pandas as pd
import numpy as np
import xgboost as xgb
import shap
from sklearn.metrics import root_mean_squared_error, r2_score
import warnings
warnings.filterwarnings('ignore')

# 1. Load the (historical) genus-level matrix
file_path = "../../data/historical/BIOTWIN_GENUS_ML_MATRIX.csv"
print("Loading Biological Matrix to Replicate Liu et al. Methodology on Duber 2024...")
df = pd.read_csv(file_path)

# 2. Isolate Duber_2024 (Paper_ID label; reported in-text as Duber et al. 2025)
study_id = 'Duber_2024'
df_local = df[df['Paper_ID'] == study_id].copy()
print(f"Isolated {study_id}. Total samples: {len(df_local)}")

# 3. Clean sparse columns only -- NOT an explicit downstream-co-product exclusion
cols_to_drop = ['succinate', 'Lactose', 'lactose', 'i-propanol', 'izo-butanol',
                 'NAOH', 'naoh', 'Paper_ID', 'Sample_ID', 'Groups']
df_local = df_local.drop(columns=[c for c in cols_to_drop if c in df_local.columns])

TARGET_METABOLITE = 'Caproate'

# 4. The Liu et al. validation split (train on Reactor B1, test on Reactor B2)
df_local['BIOREACTOR'] = df_local['BIOREACTOR'].astype(str).str.strip()
df_train = df_local[df_local['BIOREACTOR'] == 'B1']
df_test = df_local[df_local['BIOREACTOR'] == 'B2']
print(f"Training on Bioreactor B1 (Samples: {len(df_train)}) | Testing on Bioreactor B2 (Samples: {len(df_test)})")

X_train = df_train.drop(columns=[TARGET_METABOLITE, 'BIOREACTOR']).apply(pd.to_numeric, errors='coerce').dropna(axis=1, how='all')
y_train = pd.to_numeric(df_train[TARGET_METABOLITE], errors='coerce').fillna(0)

X_test = df_test.drop(columns=[TARGET_METABOLITE, 'BIOREACTOR']).apply(pd.to_numeric, errors='coerce').dropna(axis=1, how='all')
y_test = pd.to_numeric(df_test[TARGET_METABOLITE], errors='coerce').fillna(0)

X_test = X_test[X_train.columns]  # align columns; no scaling applied at this stage

xgb_params = {
    'n_estimators': 150, 'learning_rate': 0.05, 'max_depth': 3,
    'subsample': 0.8, 'colsample_bytree': 0.8, 'random_state': 42, 'n_jobs': -1,
}

print("\n--- TRAINING INTRA-STUDY MODEL (LIU ET AL. REPLICATION) ---")
model = xgb.XGBRegressor(**xgb_params).fit(X_train, y_train)

y_pred = model.predict(X_test)
r2 = r2_score(y_test, y_pred)
rmse = root_mean_squared_error(y_test, y_pred)

print("-" * 60)
print("THE SINGLE-EXPERIMENT ILLUSION RESULTS (DUBER 2024)")
print(f"Intra-Study R-squared (R2): {r2:.3f}")
print(f"Intra-Study RMSE: {rmse:.3f} mM C")
print("-" * 60)

explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_train)
mean_shap = np.abs(shap_values).mean(axis=0)
shap_df = pd.DataFrame({'Feature': X_train.columns, 'SHAP_Importance': mean_shap}).sort_values(
    by='SHAP_Importance', ascending=False)

print("\n--- TOP 10 LOCAL FEATURES ---")
print(shap_df.head(10).to_string(index=False))

# NOTE ON MANUSCRIPT TEXT: the current draft describes this run as excluding
# downstream co-products. The SHAP ranking above will show Valerate and
# Heptanoate among the top features, confirming they were not excluded here.
# Suggested correction: "R2 = -0.228 under an earlier XGBoost configuration
# lacking both explicit downstream-co-product exclusion and confirmed
# training-only scaling" in place of "XGBoost without downstream co-products."
