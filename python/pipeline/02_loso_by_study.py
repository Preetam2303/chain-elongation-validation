# 02_loso_by_study.py
# Leave-One-Study-Out cross-validation (Methods 2.6, scheme ii), run across
# four incremental feature tiers. Tier 4 (all features included) is the
# result reported as "LOSO, full genus set" in Table 3b; all four tiers
# together are Table 4 and Figure S1. This is the single script that answers
# "does feature curation rescue cross-study generalization" -- it doesn't,
# at any tier (Section 3.3).
#
# Feature lists are built via sorted(set(...)) deliberately: Python's default
# set iteration order is process-randomized, and an unsorted feature list can
# silently change which columns XGBoost's colsample_bytree samples at a given
# index even under a fixed random_state. This was a real reproducibility bug
# caught during this project -- see the project log / commit history.

import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, root_mean_squared_error
import warnings
warnings.filterwarnings('ignore')

MATRIX_PATH = "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"

df = pd.read_csv(MATRIX_PATH)

TARGET = 'Caproate'
papers = df['Paper_ID'].astype(str).str.strip()
y = pd.to_numeric(df[TARGET], errors='coerce').fillna(0)

cat_cols = ['Feed_Complexity', 'Primary_Carbon_Signature']
for c in cat_cols:
    if c in df.columns:
        df[c] = df[c].astype('category')

genus_cols = sorted([c for c in df.columns if c.startswith('g__')])
op_params = ['PH', 'TEMP', 'HRT']
foundational_feeds = ['Lactate', 'Acetate', 'Ethanol']
intermediates = ['Butyrate', 'Valerate', 'Isovalerate', 'Propionate']
downstream_proxies = ['Caprylate', 'Heptanoate']

step_1_base = op_params + cat_cols + genus_cols
step_2_feeds = step_1_base + foundational_feeds
step_3_intermediates = step_2_feeds + intermediates
step_4_downstream = step_3_intermediates + downstream_proxies

ladders = {
    "1. Base Ops + Genera (No Chems)": step_1_base,
    "2. + Base Feeds (Lactate/Acetate/Ethanol)": step_2_feeds,
    "3. + Intermediates (Butyrate/Valerate/etc.)": step_3_intermediates,
    "4. + Downstream Proxies (Caprylate/Heptanoate) -- reported as Table 3b 'LOSO, full genus set'": step_4_downstream,
}

xgb_params = {
    'n_estimators': 150, 'learning_rate': 0.05, 'max_depth': 3,
    'subsample': 0.8, 'colsample_bytree': 0.8, 'random_state': 42,
    'n_jobs': -1, 'enable_categorical': True, 'tree_method': 'hist',
}

print("--- RUNNING ZERO-LEAKAGE LOSO SUBSTRATE LADDER (Table 4 / Figure S1) ---")
print("-" * 75)

logo = LeaveOneGroupOut()

for ladder_name, feature_list in ladders.items():
    valid_features = sorted(list(set([f for f in feature_list if f in df.columns])))
    X = df[valid_features].copy()

    numeric_cols = [c for c in X.columns if c not in cat_cols]
    for c in numeric_cols:
        X[c] = pd.to_numeric(X[c], errors='coerce').fillna(0)

    r2_scores, rmse_scores = [], []

    for train_idx, test_idx in logo.split(X, y, groups=papers):
        X_train, X_test = X.iloc[train_idx].copy(), X.iloc[test_idx].copy()
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        # STRICT LEAKAGE FIX: scaler fit on training partition only
        scaler = StandardScaler()
        if len(numeric_cols) > 0:
            X_train[numeric_cols] = scaler.fit_transform(X_train[numeric_cols])
            X_test[numeric_cols] = scaler.transform(X_test[numeric_cols])

        model = xgb.XGBRegressor(**xgb_params).fit(X_train, y_train)
        y_pred = model.predict(X_test)
        if len(y_test) > 1:
            r2_scores.append(r2_score(y_test, y_pred))
            rmse_scores.append(root_mean_squared_error(y_test, y_pred))

    print(f"{ladder_name:75s}")
    print(f"    R2: {np.mean(r2_scores):6.3f} (\u00b1 {np.std(r2_scores):.3f}) | "
          f"RMSE: {np.mean(rmse_scores):6.1f} (\u00b1 {np.std(rmse_scores):.1f}) mM C")

print("-" * 75)
