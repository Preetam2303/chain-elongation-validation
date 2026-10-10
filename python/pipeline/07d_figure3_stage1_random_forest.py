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
# Published output (Duber_2024 B1 -> B2: R2 0.217) came from the historical
# matrix with rows paired by row order. Since 2026-10-10 the pairs are named
# explicitly (pipeline_settings.REACTOR_PAIRS, the same pairs as Table 5),
# unmeasured cells stay blank instead of 0 (scikit-learn's random forest
# handles blanks), and the matrix is the run's Illumina genus matrix from 00
# (BIOTWIN_GENUS), without inoculum samples. The broad predictor list is kept
# on purpose: Figure 3 shows what this earliest, loosest setup gives.
#
# (Paper_ID "Duber_2024" is the internal pipeline label for what is cited
# throughout the manuscript as Duber et al., 2025 -- see Table 1.)

import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error
import warnings

from pipeline_settings import REACTOR_PAIRS

warnings.filterwarnings('ignore')

# 1. Load the historical genus-level matrix
file_path = os.environ.get("BIOTWIN_GENUS", "../../data/historical/BIOTWIN_GENUS_ML_MATRIX.csv")
print("Loading Genus-Level Matrix for Paper-by-Paper Evaluation...")
df = pd.read_csv(file_path)

results = []

# 2. Evaluate reactor-to-reactor transferability for each named pair
df['Paper_ID'] = df['Paper_ID'].astype(str).str.strip()
df['BIOREACTOR'] = df['BIOREACTOR'].astype(str).str.strip()
cols_to_drop = ['succinate', 'Lactose', 'lactose', 'i-propanol', 'izo-butanol',
                'NAOH', 'naoh', 'Paper_ID', 'Sample_ID', 'BIOREACTOR', 'Groups',
                'Op_Mode_Binary']
for paper, (br_train, br_test) in REACTOR_PAIRS.items():
    df_p = df[df['Paper_ID'] == paper]
    if df_p.empty:
        continue  # e.g. the Nanopore study, which is not in the Illumina genus matrix
    train_df = df_p[df_p['BIOREACTOR'] == br_train]
    test_df = df_p[df_p['BIOREACTOR'] == br_test]

    X_train = train_df.drop(columns=[c for c in cols_to_drop if c in train_df.columns] + ['Caproate'],
                            errors='ignore').select_dtypes(include=[np.number])
    y_train = pd.to_numeric(train_df['Caproate'], errors='coerce')
    X_test = test_df.drop(columns=[c for c in cols_to_drop if c in test_df.columns] + ['Caproate'],
                          errors='ignore').select_dtypes(include=[np.number]).reindex(columns=X_train.columns)
    y_test = pd.to_numeric(test_df['Caproate'], errors='coerce')

    # Train and Evaluate Model -- no scaling, no explicit leakage exclusion
    model = RandomForestRegressor(n_estimators=150, random_state=42)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    results.append({
        'Paper': paper,
        'Train_Reactor': br_train,
        'Test_Reactor': br_test,
        'Train_Size': len(train_df),
        'Test_Size': len(test_df),
        'R2': r2_score(y_test, y_pred),
        'RMSD': np.sqrt(mean_squared_error(y_test, y_pred)),
    })

res_df = pd.DataFrame(results)
print("\n" + "=" * 80)
print("PAPER-BY-PAPER INTRA-STUDY GENERALIZATION RESULTS")
print("=" * 80)
print(res_df.to_string(index=False))
print("=" * 80)
