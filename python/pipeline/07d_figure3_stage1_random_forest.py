# 07d_figure3_stage1_random_forest.py
# Figure 3, stage 1: the earliest, least leakage-controlled evaluation of the
# same reactor-pair transfer, using a plain sklearn RandomForestRegressor with
# no scaling and no explicit downstream-co-product exclusion (only a short
# list of sparse/junk columns is dropped). This is the starting point of the
# three-stage progression reported in Figure 3 and Section 3.5: R2 = 0.217 ->
# -0.228 (07c) -> -5.869 (05), as leakage controls were progressively
# tightened on the identical Duber et al. (2025) B1->B2 samples.
#
# Requires the historical BIOTWIN_GENUS_ML_MATRIX.csv snapshot (see
# ../../data/historical/), the same file 07c depends on -- this and 07c are
# the two scripts in this repository that are not reproducible from
# 01_genus_aggregation.py's current output.
#
# Expected output (confirmed against the original run, and independently
# cross-checked against every prior record of these four values throughout
# this project's development):
#                          Paper Train_Reactor Test_Reactor  Train_Size  Test_Size       R2      RMSD
#                 Brodowski_2025            B1           B2           7          7 0.455090 55.564719
#                     Duber_2024            B1           B2          30         21 0.217206 50.892162
#                     Duber_2022            R1           R2           3          3 0.003083 43.651106
# Brodowski_2022_EXTERNAL_ACETATE            B1           B2           9          7 0.506065 91.911224
#
# (Paper_ID "Duber_2024" is the internal pipeline label for what is cited
# throughout the manuscript as Duber et al., 2025 -- see Table 1.)

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error
import warnings
warnings.filterwarnings('ignore')

# 1. Load the historical genus-level matrix
file_path = "../../data/historical/BIOTWIN_GENUS_ML_MATRIX.csv"
print("Loading Genus-Level Matrix for Paper-by-Paper Evaluation...")
df = pd.read_csv(file_path)

results = []

# 2. Loop through each paper and evaluate reactor-to-reactor transferability
for paper in df['Paper_ID'].unique():
    df_p = df[df['Paper_ID'] == paper].copy()
    bioreactors = df_p['BIOREACTOR'].dropna().unique()

    # Check if the paper has at least 2 distinct bioreactors/systems
    if len(bioreactors) >= 2:
        br_train = bioreactors[0]
        br_test = bioreactors[1]

        train_df = df_p[df_p['BIOREACTOR'].astype(str).str.strip() == str(br_train).strip()]
        test_df = df_p[df_p['BIOREACTOR'].astype(str).str.strip() == str(br_test).strip()]

        if len(train_df) >= 3 and len(test_df) >= 3:
            cols_to_drop = ['succinate', 'Lactose', 'lactose', 'i-propanol', 'izo-butanol',
                             'NAOH', 'naoh', 'Paper_ID', 'Sample_ID', 'BIOREACTOR', 'Groups',
                             'Op_Mode_Binary']

            X_train = train_df.drop(columns=[c for c in cols_to_drop if c in train_df.columns] + ['Caproate'],
                                     errors='ignore').select_dtypes(include=[np.number]).fillna(0)
            y_train = pd.to_numeric(train_df['Caproate'], errors='coerce').fillna(0)

            X_test = test_df.drop(columns=[c for c in cols_to_drop if c in test_df.columns] + ['Caproate'],
                                   errors='ignore').select_dtypes(include=[np.number]).fillna(0)
            y_test = pd.to_numeric(test_df['Caproate'], errors='coerce').fillna(0)

            # Align features safely
            X_test = X_test.reindex(columns=X_train.columns, fill_value=0)

            # Train and Evaluate Model -- no scaling, no explicit leakage exclusion
            model = RandomForestRegressor(n_estimators=150, random_state=42)
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)

            r2 = r2_score(y_test, y_pred)
            rmsd = np.sqrt(mean_squared_error(y_test, y_pred))

            results.append({
                'Paper': paper,
                'Train_Reactor': br_train,
                'Test_Reactor': br_test,
                'Train_Size': len(train_df),
                'Test_Size': len(test_df),
                'R2': r2,
                'RMSD': rmsd,
            })

res_df = pd.DataFrame(results)
print("\n" + "=" * 80)
print("PAPER-BY-PAPER INTRA-STUDY GENERALIZATION RESULTS")
print("=" * 80)
print(res_df.to_string(index=False))
print("=" * 80)
