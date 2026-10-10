# 07b_algorithm_comparison_lightgbm_standard.py
# The second half of the algorithm-selection comparison in Methods 2.7: a
# standard gradient-boosted (gbdt) LightGBM configuration against XGBoost on
# the exact same feature set and 80/20 split (random_state=42) as the Table 3a
# single-split baseline (10_table3a_naive_splits.py), for direct comparison.
#
# Feature set: the core set (pH, temperature, HRT, lactate, acetate, ethanol,
# butyrate, propionate, the 2 categoricals) + ASVs, as everywhere else
# (pipeline_settings.py). Both models use the shared settings (150 trees,
# learning rate 0.05, depth 3, no subsampling), so the comparison is like for
# like. Unmeasured cells stay blank; both libraries handle blanks.
#
# Published before this change: XGBoost 0.915 (typed in, not computed),
# LightGBM 0.861, from a feature set that included NaOH and DAY.

import os
import warnings

import pandas as pd
import xgboost as xgb
from lightgbm import LGBMRegressor
from sklearn.metrics import r2_score, root_mean_squared_error
from sklearn.model_selection import train_test_split

from pipeline_settings import CAT_COLS, CORE_CHEM, OPS, TARGET, XGB_PARAMS

warnings.filterwarnings('ignore')

# ==============================================================================
# 1. LOAD THE PRUNED MATRIX & ENVIRONMENT SETUP
# ==============================================================================
file_path = os.environ.get("BIOTWIN_PRUNED", "../../data/historical/BIOTWIN_PRUNED_ML_MATRIX.csv")
print("Loading Pruned ML Matrix...")
df = pd.read_csv(file_path, low_memory=False)

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

X = df[core_cols + core_asvs]
y = df[TARGET]
# The same 80/20 split as the Table 3a single-split baseline
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

print(f"Training Matrix built with {X_train.shape[1]} features.")
print(f"Training on {X_train.shape[0]} samples, testing on {X_test.shape[0]} samples.")

# ==============================================================================
# 3. XGBOOST (shared settings) AND STANDARD GRADIENT-BOOSTED LIGHTGBM
# ==============================================================================
xgb_pred = xgb.XGBRegressor(**XGB_PARAMS).fit(X_train, y_train).predict(X_test)

# 'gbdt' (Gradient Boosting Decision Tree), not 'rf' -- this is the standard
# boosting comparison, distinct from 07a's bagged Random-Forest-mode run.
lgbm_ultimate = LGBMRegressor(
    boosting_type='gbdt', n_estimators=150, learning_rate=0.05, max_depth=3,
    num_leaves=7,             # capped safely below 2^max_depth
    min_child_samples=5,      # minimum samples in a leaf to allow a split
    random_state=42, n_jobs=-1, verbose=-1,
)
print("\nTraining Standard LightGBM (GBDT) on the core set + ASVs...")
lgbm_pred = lgbm_ultimate.fit(X_train, y_train).predict(X_test)

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
