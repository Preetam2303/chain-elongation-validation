# 07b_algorithm_comparison_lightgbm_standard.py
# The second half of the algorithm-selection comparison in Methods 2.7: a
# standard gradient-boosted (gbdt) LightGBM configuration on the exact same
# Ultimate feature subset and 80/20 split (random_state=42) as the XGBoost
# single-split baseline (Table 3a), for direct comparison.
#
# Expected output (confirmed against the original run):
#   XGBoost baseline R2:  0.809
#   LightGBM (gbdt) R2:   0.777, RMSE = 52.696 mM C
#
# XGBoost's level-wise tree growth outperformed both this and the
# Random-Forest-mode configuration (07a) under identical conditions and was
# adopted throughout the rest of this pipeline.

import pandas as pd
import numpy as np
from lightgbm import LGBMRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import root_mean_squared_error, r2_score

# ==============================================================================
# 1. LOAD THE PRUNED MATRIX & ENVIRONMENT SETUP
# ==============================================================================
file_path = "../../data/historical/BIOTWIN_PRUNED_ML_MATRIX.csv"
print("Loading Pruned ML Matrix...")
df = pd.read_csv(file_path)

core_asvs = [col for col in df.columns if str(col).startswith('ASV_')]

# ==============================================================================
# 2. DEFINE THE ULTIMATE FEATURE SUBSET & TARGET
# ==============================================================================
TARGET_METABOLITE = 'Caproate'
op_cols = ['PH', 'TEMP', 'HRT', 'NAOH', 'DAY']
substrate_cols = ['Lactate', 'Acetate', 'Butyrate', 'Ethanol', 'Propionate']
ultimate_features = op_cols + substrate_cols + core_asvs

# Fix Target NaNs ONLY (Day 0 samples where Caproate is 0)
df[TARGET_METABOLITE] = pd.to_numeric(df[TARGET_METABOLITE], errors='coerce').fillna(0)
for col in op_cols + substrate_cols:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# ==============================================================================
# 3. CONSTRUCT MATRIX AND TRAIN/TEST SPLIT
# ==============================================================================
X = df[ultimate_features]
y = df[TARGET_METABOLITE]
# Replicates the exact 80/20 train/test split used for the XGBoost baseline
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

print(f"Training Matrix built with {X_train.shape[1]} features.")
print(f"Training on {X_train.shape[0]} samples, testing on {X_test.shape[0]} samples.")

# ==============================================================================
# 4. INITIALIZE STANDARD GRADIENT BOOSTED LIGHTGBM
# ==============================================================================
# 'gbdt' (Gradient Boosting Decision Tree), not 'rf' -- this is the standard
# boosting comparison, distinct from 07a's bagged Random-Forest-mode run.
lgbm_ultimate = LGBMRegressor(
    boosting_type='gbdt', n_estimators=150, learning_rate=0.05, max_depth=3,
    num_leaves=7,             # capped safely below 2^max_depth
    min_child_samples=5,      # minimum samples in a leaf to allow a split
    subsample=0.8, colsample_bytree=0.8,
    random_state=42, n_jobs=-1, verbose=-1,
)

# ==============================================================================
# 5. TRAIN AND EVALUATE
# ==============================================================================
print("\nTraining Standard LightGBM (GBDT) on Ultimate Feature Subset...")
lgbm_ultimate.fit(X_train, y_train)
y_pred = lgbm_ultimate.predict(X_test)

r2 = r2_score(y_test, y_pred)
rmse = root_mean_squared_error(y_test, y_pred)

print("\n" + "=" * 50)
print("STANDARD LIGHTGBM ULTIMATE MODEL RESULTS")
print("=" * 50)
print(f"Target Variable:     {TARGET_METABOLITE}")
print(f"XGBoost Baseline R2: 0.809")
print(f"LightGBM R-squared:  {r2:.3f}")
print(f"LightGBM RMSE:       {rmse:.3f} mM C")
print("=" * 50)
