# 07a_algorithm_comparison_lightgbm_rf.py
# Algorithm-selection comparison referenced in Methods 2.7: tests a bagged,
# Random-Forest-mode LightGBM configuration against three feature scopes,
# on the same ASV-level 80/20 single-split baseline used elsewhere in early
# algorithm selection (random_state=42, identical split across all algorithm
# comparisons in this project). This is the run that surfaced the concrete
# target-leakage demonstration cited in Methods 2.7: admitting the full
# chemistry panel, including the downstream co-products caprylate, heptanoate
# and isovalerate, nearly doubled apparent accuracy (0.359 -> 0.657 in the
# published run) on an otherwise identical feature scope, motivating their
# exclusion from every leakage-controlled script elsewhere in this pipeline
# (Methods 2.5).
#
# Feature scopes follow the decision of 2026-10-10 (pipeline_settings.py):
#   1. Base:     pH, temperature, HRT, the 2 categoricals + ASVs
#   2. Ultimate: the core set (adds lactate, acetate, ethanol, butyrate,
#                propionate) + ASVs
#   3. Maximum:  the core set + the Illumina-only compounds (valerate,
#                isovalerate, caprylate, heptanoate, isocaproate) + ASVs, the
#                leakage demonstration
# NaOH, DAY and the other never-modelled columns are no longer used; blank
# (unmeasured) cells stay blank. The values published before this change
# (0.354 / 0.359 / 0.657) came from the earlier scopes, which included NaOH and
# DAY. LightGBM's Random-Forest mode needs bagging (subsample < 1), so this
# script keeps it; it is the only model in the pipeline that samples rows.
#
# Note: this operates on the ASV-level (pre-genus-collapse) matrix, not the
# 184-sample genus-level matrix used throughout the rest of this pipeline --
# it predates the genus aggregation step and is retained here exactly as run,
# for methodological traceability rather than as a headline result.

import os
import pandas as pd
import numpy as np
from lightgbm import LGBMRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import root_mean_squared_error, r2_score
import shap
import copy

from pipeline_settings import CAT_COLS, CORE_CHEM, ILLUMINA_ONLY_TIER3, ILLUMINA_ONLY_TIER4, OPS, TARGET

# ==============================================================================
# 1. LOAD THE PRUNED MATRIX & IDENTIFY ASVs
# ==============================================================================
file_path = os.environ.get("BIOTWIN_PRUNED", "../../data/historical/BIOTWIN_PRUNED_ML_MATRIX.csv")
print("Loading Pruned ML Matrix...")
df = pd.read_csv(file_path)

core_asvs = [col for col in df.columns if str(col).startswith('ASV_')]
print(f"Loaded {df.shape[0]} samples and {len(core_asvs)} core ASVs.")

# ==============================================================================
# 2. PREPROCESS: numeric predictors, unmeasured cells left blank
# ==============================================================================
TARGET_METABOLITE = TARGET
base_cols = OPS + CAT_COLS
core_cols = OPS + CORE_CHEM + CAT_COLS
maximum_cols = core_cols + ILLUMINA_ONLY_TIER3 + ILLUMINA_ONLY_TIER4

for col in OPS + CORE_CHEM + ILLUMINA_ONLY_TIER3 + ILLUMINA_ONLY_TIER4 + [TARGET]:
    df[col] = pd.to_numeric(df[col], errors='coerce')
for col in CAT_COLS:
    df[col] = df[col].astype('category')
if df[TARGET].isna().any():
    raise SystemExit("Caproate is missing in the ASV-level matrix; every modelled row needs it.")

# ==============================================================================
# 3. SETUP THE LIGHTGBM RANDOM FOREST
# ==============================================================================
models_to_test = {
    "1. Base Model (Ops + categoricals + ASVs)": base_cols + core_asvs,
    "2. Ultimate Model (core set + ASVs)": core_cols + core_asvs,
    "3. Maximum Model (core set + Illumina-only compounds + ASVs)": maximum_cols + core_asvs,
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
mean_abs = pd.Series(np.abs(shap_values.values).mean(axis=0), index=best_X_test.columns)
# Only the order is printed: for LightGBM's Random-Forest mode the SHAP values
# are summed over trees rather than averaged, so their size is not in mM C.
print(f"Top 10 features by mean |SHAP| ({best_name}):")
for rank, feature in enumerate(mean_abs.sort_values(ascending=False, kind='mergesort').head(10).index, 1):
    print(f"  {rank:2d}. {feature}")
print("If the Maximum model wins and the Illumina-only compounds (caprylate, heptanoate, "
      "isovalerate) lead, that is the target-leakage demonstration cited in Methods 2.7.")
