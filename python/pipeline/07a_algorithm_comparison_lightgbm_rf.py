# 07a_algorithm_comparison_lightgbm_rf.py
# Algorithm-selection comparison referenced in Methods 2.7: tests a bagged,
# Random-Forest-mode LightGBM configuration against three feature scopes,
# on the same ASV-level 80/20 single-split baseline used elsewhere in early
# algorithm selection (random_state=42, identical split across all algorithm
# comparisons in this project). This is the run that surfaced the concrete
# target-leakage demonstration cited in Methods 2.7: the winning "Maximum"
# feature scope relies most heavily on isovalerate, heptanoate, and caprylate
# -- downstream co-products of caproate -- as its top predictors by SHAP
# attribution, motivating their exclusion from every leakage-controlled
# script elsewhere in this pipeline (Methods 2.5).
#
# Expected output (confirmed against the original run):
#   1. Base Model (Ops + ASVs):                    R2 = 0.355, RMSE = 89.550 mM C
#   2. Ultimate Model (Ops + Limited Subs + ASVs):  R2 = 0.381, RMSE = 87.781 mM C
#   3. Maximum Model (Ops + All Chems + ASVs):      R2 = 0.415, RMSE = 85.312 mM C
#   Winner: Maximum Model, R2 = 0.415 -- still well below XGBoost's 0.809 on
#   the equivalent feature scope (see 02_loso_by_study.py / Table 3a).
#
# Note: this operates on the ASV-level (pre-genus-collapse) matrix, not the
# 184-sample genus-level matrix used throughout the rest of this pipeline --
# it predates the genus aggregation step and is retained here exactly as run,
# for methodological traceability rather than as a headline result.

import pandas as pd
import numpy as np
from lightgbm import LGBMRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import root_mean_squared_error, r2_score
import shap
import copy

# ==============================================================================
# 1. LOAD THE PRUNED MATRIX & IDENTIFY ASVs
# ==============================================================================
file_path = "../../data/historical/BIOTWIN_PRUNED_ML_MATRIX.csv"
print("Loading Pruned ML Matrix...")
df = pd.read_csv(file_path)

core_asvs = [col for col in df.columns if str(col).startswith('ASV_')]
print(f"Loaded {df.shape[0]} samples and {len(core_asvs)} core ASVs.")

# ==============================================================================
# 2. PREPROCESS (FIX ONLY CAPROATE NaNs)
# ==============================================================================
TARGET_METABOLITE = 'Caproate'
op_cols = ['PH', 'TEMP', 'HRT', 'NAOH', 'DAY']
substrate_cols = ['Lactate', 'Acetate', 'Butyrate', 'Ethanol', 'Propionate']
all_chem_cols = ['Lactate', 'Lactose', 'succinate', 'Acetate', 'Propionate',
                  'Isobutyrate', 'Butyrate', 'Isovalerate', 'Valerate',
                  'Isocaproate', 'Heptanoate', 'Caprylate', 'Ethanol',
                  'i-propanol', 'Propanol', 'izo-butanol', 'Butanol']

df[TARGET_METABOLITE] = pd.to_numeric(df[TARGET_METABOLITE], errors='coerce').fillna(0)
for col in op_cols + all_chem_cols:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# ==============================================================================
# 3. SETUP THE LIGHTGBM RANDOM FOREST
# ==============================================================================
models_to_test = {
    "1. Base Model (Ops + ASVs)": op_cols + core_asvs,
    "2. Ultimate Model (Ops + Limited Subs + ASVs)": op_cols + substrate_cols + core_asvs,
    "3. Maximum Model (Ops + All Chems + ASVs)": op_cols + all_chem_cols + core_asvs,
}

# 'rf' boosting_type requires both bagging and feature sampling to be < 1.0
lgbm_rf_base = LGBMRegressor(
    boosting_type='rf', n_estimators=150, max_depth=5,
    subsample=0.8, subsample_freq=1, colsample_bytree=0.8,
    random_state=42, n_jobs=-1, verbose=-1,
)

best_model, best_name, best_r2, best_X_test = None, "", -float('inf'), None

# ==============================================================================
# 4. TRAIN AND EVALUATE THE 3 MODELS
# ==============================================================================
print("\n=== LIGHTGBM RANDOM FOREST COMPARISON ===")
for model_name, feature_list in models_to_test.items():
    X = df[feature_list]
    y = df[TARGET_METABOLITE]
    # random_state=42 matches the XGBoost single-split baseline (Table 3a)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    rf_model = copy.deepcopy(lgbm_rf_base)
    rf_model.fit(X_train, y_train)
    y_pred = rf_model.predict(X_test)

    r2 = r2_score(y_test, y_pred)
    rmse = root_mean_squared_error(y_test, y_pred)

    print(f"\n{model_name}:")
    print(f"  Features: {len(feature_list)}")
    print(f"  R-squared (R2): {r2:.3f}")
    print(f"  RMSE: {rmse:.3f} mM C")

    if r2 > best_r2:
        best_r2, best_name, best_model, best_X_test = r2, model_name, copy.deepcopy(rf_model), X_test.copy()

print("\n" + "=" * 50)
print(f"WINNING RF MODEL: {best_name}")
print(f"BEST R2 SCORE:    {best_r2:.3f}")
print("=" * 50)

# ==============================================================================
# 5. SHAP FOR THE WINNER -- surfaces the downstream-co-product leakage pattern
# ==============================================================================
print("\nCalculating SHAP values for the winning model...")
explainer = shap.TreeExplainer(best_model)
shap_values = explainer(best_X_test)
print("SHAP values calculated. Top features by mean |SHAP| for the Maximum model "
      "should show isovalerate, heptanoate, and caprylate leading -- the target-"
      "leakage demonstration cited in Methods 2.7.")
