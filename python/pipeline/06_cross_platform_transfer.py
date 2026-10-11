# 06_cross_platform_transfer.py
# Table 3b, platform transfer: train on all Illumina rows, test on the Nanopore
# study, with the core predictor set and the settings in pipeline_settings.py
# (no subsampling). Run twice: genera as raw read counts and as relative
# abundance. The training-mean benchmark (predict the Illumina rows' mean
# caproate) is printed next to them (review item 5).
#
# Spearman rho between measured and predicted caproate is computed on
# sampling-day means, so the Nanopore triplicates count once (18 days).

import os
import warnings

import xgboost as xgb

from pipeline_settings import (CAT_COLS, CORE_CHEM, NANOPORE_ID, OPS, TARGET, XGB_PARAMS, load_matrix,
                               present, r2_rmse, sampling_units, save_table, spearman_by_unit,
                               training_mean_r2_rmse)

warnings.filterwarnings('ignore')

file_path = os.environ.get("BIOTWIN_MATRIX", "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")

table = []
for relative, label in [(False, "raw counts"), (True, "relative abundance")]:
    df, genus_cols = load_matrix(file_path, relative=relative)
    units = sampling_units(df)
    train_df = df[df['Paper_ID'] != NANOPORE_ID]
    test_df = df[df['Paper_ID'] == NANOPORE_ID]
    features = sorted(set(genus_cols + present(OPS + CORE_CHEM + CAT_COLS, df)))

    model = xgb.XGBRegressor(**XGB_PARAMS).fit(train_df[features], train_df[TARGET])
    y_pred = model.predict(test_df[features])
    r2, rmse = r2_rmse(test_df[TARGET], y_pred)
    rho, n_days = spearman_by_unit(test_df[TARGET], y_pred, units[test_df.index])

    print(f"--- {label.upper()}: ILLUMINA -> NANOPORE TRANSFER ---")
    print(f"Training Samples (Illumina): {len(train_df)}")
    print(f"Testing Samples (Nanopore): {len(test_df)} ({n_days} sampling days)")
    print(f"Features: {len(features)}")
    print(f"Transfer R2:   {r2:.3f}")
    print(f"Transfer RMSE: {rmse:.3f} mM C")
    print(f"Spearman rho (sampling-day means): {rho:.3f}")
    table.append({'table': 'Table 3b', 'scheme': f'Platform transfer, {label}', 'model': 'XGBoost, core set',
                  'train_n': len(train_df), 'test_n': len(test_df), 'r2': r2, 'rmse': rmse,
                  'spearman_rho_day_means': rho})

b_r2, b_rmse = training_mean_r2_rmse(train_df[TARGET], test_df[TARGET])
print("--- TRAINING-MEAN BENCHMARK (predict the Illumina mean) ---")
print(f"Benchmark R2:   {b_r2:.3f}")
print(f"Benchmark RMSE: {b_rmse:.3f} mM C")
table.append({'table': 'Table 3b', 'scheme': 'Platform transfer', 'model': 'Training-mean benchmark',
              'train_n': len(train_df), 'test_n': len(test_df), 'r2': b_r2, 'rmse': b_rmse,
              'spearman_rho_day_means': float('nan')})
save_table(table, "06_platform_transfer.csv")
