# 05_intrastudy_transferability_table5.py
# Table 5 (and Figure 3, stage 3): reactor-to-reactor transfer within a study.
# Train on one reactor, test on the other reactor of the same study, with the
# core predictor set, training-reactor-only scaling and no subsampling
# (pipeline_settings.py), so every number reproduces exactly.
#
# The pairs are named explicitly (pipeline_settings.REACTOR_PAIRS; review
# item 9), not taken from row order. Studies without two reactors of at least
# 3 samples each are listed with the reason (NOT_IN_TABLE5).
#
# The two reactors of each pair did NOT run under matched conditions: the
# script prints each reactor's feed descriptors and its pH, temperature and
# HRT values, so the table can say how they differed.
#
# Spearman rho between measured and predicted caproate is computed on
# sampling-unit means: the three same-day replicate samples of Duber 2024 and
# the Nanopore study count once (decision 2026-10-07).

import os
import warnings

import numpy as np
import pandas as pd
import xgboost as xgb

from pipeline_settings import (CAT_COLS, CORE_CHEM, NOT_IN_TABLE5, OPS, REACTOR_PAIRS, TARGET, XGB_PARAMS,
                               load_matrix, present, r2_rmse, sampling_units, save_table, scale,
                               spearman_by_unit)

warnings.filterwarnings('ignore')

file_path = os.environ.get("BIOTWIN_MATRIX", "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")
print("Loading Grand Merge Matrix for Clean Intra-Study Evaluation...")
df, genus_cols = load_matrix(file_path)
units = sampling_units(df)
valid_features = sorted(set(present(OPS + CORE_CHEM + CAT_COLS, df) + genus_cols))


def counts(values):
    v = pd.Series(values).dropna()
    if v.empty:
        return "blank"
    vc = v.value_counts().sort_index()
    return ", ".join(f"{k:g} (x{n})" if isinstance(k, (int, float, np.floating)) else f"{k} (x{n})"
                     for k, n in vc.items())


print(f"\n--- REACTOR CONDITIONS PER PAIR (the two reactors ran differently by design) ---")
for paper, (br_train, br_test) in REACTOR_PAIRS.items():
    print(f"\n{paper}")
    for br in (br_train, br_test):
        d = df[(df['Paper_ID'] == paper) & (df['BIOREACTOR'] == br)]
        feed = ", ".join(sorted({f"{a} / {b}" for a, b in zip(d['Feed_Complexity'], d['Primary_Carbon_Signature'])}))
        print(f"  {br} ({len(d)} samples) | feed: {feed}")
        for c in OPS:
            print(f"      {c:<5s} {counts(d[c])}")

results = []
print(f"\n--- EXECUTING LEAK-PROOF INTRA-STUDY GENERALIZATION (XGBOOST) ---")
print("-" * 88)
for paper, (br_train, br_test) in REACTOR_PAIRS.items():
    train = df[(df['Paper_ID'] == paper) & (df['BIOREACTOR'] == br_train)]
    test = df[(df['Paper_ID'] == paper) & (df['BIOREACTOR'] == br_test)]
    if len(train) < 3 or len(test) < 3:
        raise SystemExit(f"{paper}: {br_train} has {len(train)} and {br_test} has {len(test)} samples; need 3 or more each.")

    X_train, X_test = scale(train[valid_features], test[valid_features])
    model = xgb.XGBRegressor(**XGB_PARAMS).fit(X_train, train[TARGET])
    y_pred = model.predict(X_test)

    r2, rmse = r2_rmse(test[TARGET], y_pred)
    rho, n_units = spearman_by_unit(test[TARGET], y_pred, units[test.index])
    results.append({
        'Paper': paper, 'Train_Reactor': br_train, 'Test_Reactor': br_test,
        'Train_n': len(train), 'Test_n': len(test), 'Test_sampling_units': n_units,
        'R2': r2, 'RMSE (mM C)': rmse, 'Spearman_rho_unit_means': rho,
    })

res_df = pd.DataFrame(results)
print("\n" + "=" * 88)
print("CLEAN INTRA-STUDY GENERALIZATION RESULTS (NO TARGET LEAKAGE)")
print("=" * 88)
print(res_df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
print("=" * 88)
print("Not in Table 5:")
for paper, why in NOT_IN_TABLE5.items():
    if paper in set(df['Paper_ID']):
        print(f"  {paper}: {why}")
save_table(res_df, "05_table5_reactor_pairs.csv")
