import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, root_mean_squared_error
import warnings
warnings.filterwarnings('ignore')

file_path = "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"
print("Loading Grand Merge Matrix for Clean Intra-Study Evaluation...")
df = pd.read_csv(file_path)

genus_cols = sorted([col for col in df.columns if col.startswith('g__')])
df[genus_cols] = df[genus_cols].fillna(0)
row_sums = df[genus_cols].sum(axis=1).replace(0, 1)
df[genus_cols] = df[genus_cols].div(row_sums, axis=0)

TARGET = 'Caproate'
cat_cols = ['Feed_Complexity', 'Primary_Carbon_Signature']
for c in cat_cols:
    if c in df.columns:
        df[c] = df[c].astype('category')

leaking_coproducts = ['Caprylate', 'Heptanoate', 'Valerate', 'Isovalerate', 'succinate', 'Lactose', 'lactose', 'i-propanol', 'izo-butanol', 'NAOH', 'naoh']
foundational_chems = ['Lactate', 'Acetate', 'Butyrate', 'Ethanol', 'Propionate']
op_params = ['PH', 'TEMP', 'HRT']

valid_chems_ops = [c for c in foundational_chems + op_params if c in df.columns]
valid_cats = [c for c in cat_cols if c in df.columns]

valid_features = sorted(list(set(valid_chems_ops + valid_cats + genus_cols)))

xgb_params = {
    'n_estimators': 150, 'learning_rate': 0.05, 'max_depth': 3,
    'subsample': 0.8, 'colsample_bytree': 0.8, 'random_state': 42,
    'n_jobs': -1, 'enable_categorical': True, 'tree_method': 'hist'
}

results = []

print(f"\n--- EXECUTING LEAK-PROOF INTRA-STUDY GENERALIZATION (XGBOOST) ---")
print("-" * 88)

for paper in df['Paper_ID'].unique():
    df_p = df[df['Paper_ID'] == paper].copy()
    bioreactors = df_p['BIOREACTOR'].dropna().unique()
    print(f"DEBUG {paper}: bioreactors order = {list(bioreactors)}")

    if len(bioreactors) >= 2:
        br_train = bioreactors[0]
        br_test = bioreactors[1]

        train_df = df_p[df_p['BIOREACTOR'].astype(str).str.strip() == str(br_train).strip()].copy()
        test_df = df_p[df_p['BIOREACTOR'].astype(str).str.strip() == str(br_test).strip()].copy()

        if len(train_df) >= 3 and len(test_df) >= 3:
            X_train = train_df[valid_features].copy()
            y_train = pd.to_numeric(train_df[TARGET], errors='coerce').fillna(0)

            X_test = test_df[valid_features].copy()
            y_test = pd.to_numeric(test_df[TARGET], errors='coerce').fillna(0)

            numeric_cols = [c for c in valid_features if c not in valid_cats]
            for c in numeric_cols:
                X_train[c] = pd.to_numeric(X_train[c], errors='coerce').fillna(0)
                X_test[c] = pd.to_numeric(X_test[c], errors='coerce').fillna(0)

            scaler = StandardScaler()
            if len(numeric_cols) > 0:
                X_train[numeric_cols] = scaler.fit_transform(X_train[numeric_cols])
                X_test[numeric_cols] = scaler.transform(X_test[numeric_cols])

            model = xgb.XGBRegressor(**xgb_params).fit(X_train, y_train)
            y_pred = model.predict(X_test)

            r2 = r2_score(y_test, y_pred)
            rmse = root_mean_squared_error(y_test, y_pred)

            results.append({
                'Paper': paper,
                'Train_Reactor': br_train,
                'Test_Reactor': br_test,
                'Train_n': len(train_df),
                'Test_n': len(test_df),
                'R2': r2,
                'RMSE (mM C)': rmse
            })

res_df = pd.DataFrame(results)
print("\n" + "="*88)
print("CLEAN INTRA-STUDY GENERALIZATION RESULTS (NO TARGET LEAKAGE)")
print("="*88)
print(res_df.to_string(index=False, float_format=lambda x: f"{x:.3f}" if isinstance(x, float) else str(x)))
print("="*88)
