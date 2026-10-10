# pipeline_settings.py
# Settings shared by the modelling scripts (02 to 10), so they cannot drift
# apart. Decided by Preetam on 2026-10-10 (project file
# rerun/2026-10-10_pr2_plan/DECISION_headline_rows_features_2026-10-10.md):
#
#   - Core predictors in every validation scheme: pH, temperature, HRT,
#     lactate, acetate, ethanol, butyrate, propionate, the two categorical
#     descriptors, and the genera. (Figure 3's first two stages, 07d and 07c,
#     keep their broader lists on purpose.)
#   - Headline: full-matrix LOSO with exactly the core set (Table 4 full-matrix
#     tier 3 = Table 3b's LOSO row).
#   - Compounds measured only in the Illumina studies (valerate, isovalerate,
#     caprylate, heptanoate, isocaproate) enter only Table 4's Illumina-only
#     column.
#   - Never modelled: isobutyrate, propanol, butanol, i-propanol, iso-butanol,
#     succinate, NaOH, lactose. DAY is not a predictor.
#   - A compound that was not measured stays blank (NaN), never 0. Genera are
#     the exception: an absent genus really has 0 reads.
#   - One set of XGBoost settings everywhere, with no subsampling, so every
#     number reproduces exactly on any machine.
#   - Only predictors are scaled (train rows only); the target never is.
#
# Inoculum samples are never modelled (data/corrections/descriptive_only_samples.csv).
# Script 00 leaves them out of the corrected-run matrices; drop_unmodelled()
# also drops them, and the excluded samples, from any other matrix, such as
# the published ones in data/.

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import r2_score, root_mean_squared_error
from sklearn.preprocessing import StandardScaler

TARGET = 'Caproate'
NANOPORE_ID = 'Hanna_2025'

CAT_COLS = ['Feed_Complexity', 'Primary_Carbon_Signature']
OPS = ['PH', 'TEMP', 'HRT']
FEEDS = ['Lactate', 'Acetate', 'Ethanol']
INTERMEDIATES = ['Butyrate', 'Propionate']
CORE_CHEM = FEEDS + INTERMEDIATES
CORE = OPS + CORE_CHEM + CAT_COLS  # plus the genera (or ASVs)
ILLUMINA_ONLY_TIER3 = ['Valerate', 'Isovalerate']
ILLUMINA_ONLY_TIER4 = ['Caprylate', 'Heptanoate', 'Isocaproate']
NEVER_MODELLED = ['Isobutyrate', 'Propanol', 'Butanol', 'i-propanol', 'izo-butanol',
                  'succinate', 'NAOH', 'Lactose']
ID_COLS = ['Paper_ID', 'Sample_ID', 'BIOREACTOR', 'Operation_Mode', 'DAY']

XGB_PARAMS = dict(n_estimators=150, learning_rate=0.05, max_depth=3, random_state=42,
                  n_jobs=-1, enable_categorical=True, tree_method='hist')

# Same-day replicate samples (decision 2026-10-07): three samples per reactor
# and day, all paired with that day's single chemistry measurement. They are
# one sampling unit wherever rows would otherwise count as independent.
REPLICATE_STUDIES = ['Duber_2024', NANOPORE_ID]

# Table 5 reactor pairs, named explicitly (review item 9). Train on the first
# reactor, test on the second.
REACTOR_PAIRS = {
    'Brodowski_2025': ('B1', 'B2'),
    'Duber_2024': ('B1', 'B2'),
    'Brodowski_2022_EXTERNAL_ACETATE': ('B1', 'B2'),
    NANOPORE_ID: ('B1', 'B2'),
}
NOT_IN_TABLE5 = {
    'Duber_2020': 'one reactor',
    'Duber_2022': 'batch bottles, 2 per condition; no reactor with 3 or more samples',
    'Brodowski_2022': 'CSTR B1 has 1 sample (B2 has 4); the batch bottles are 1 per condition',
}


CORRECTIONS_DIR = Path(__file__).resolve().parents[2] / 'data' / 'corrections'


def drop_unmodelled(df):
    """Leave out the samples that are never modelled: the inoculum samples
    (descriptive_only_samples.csv) and the excluded ones (excluded_samples.csv).
    A no-op on the corrected-run matrices, which no longer contain them."""
    ids = set()
    for name in ('descriptive_only_samples.csv', 'excluded_samples.csv'):
        ids |= set(pd.read_csv(CORRECTIONS_DIR / name)['Sample_ID'].astype(str).str.strip())
    drop = df['Sample_ID'].astype(str).str.strip().isin(ids)
    if drop.any():
        print(f"Rows left out (inoculum or excluded samples, never modelled): {int(drop.sum())}")
    return df[~drop].reset_index(drop=True)


def load_matrix(path, relative=True):
    """Read a modelling matrix. Genera as relative abundance (default), other
    predictors numeric with unmeasured cells left blank. Returns (df, genus_cols)."""
    df = pd.read_csv(path, low_memory=False)
    for c in ['Paper_ID', 'Sample_ID', 'BIOREACTOR']:
        df[c] = df[c].astype(str).str.strip()
    df = drop_unmodelled(df)
    genus_cols = sorted(c for c in df.columns if c.startswith('g__'))
    df[genus_cols] = df[genus_cols].fillna(0)
    if relative:
        df[genus_cols] = df[genus_cols].div(df[genus_cols].sum(axis=1).replace(0, 1), axis=0)
    for c in CAT_COLS:
        if c in df.columns:
            df[c] = df[c].astype('category')
    for c in df.columns:
        if c not in genus_cols + CAT_COLS + ID_COLS:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    if df[TARGET].isna().any():
        missing = df.loc[df[TARGET].isna(), 'Sample_ID'].tolist()
        raise SystemExit(f"Caproate is missing for {missing}. Every modelled row needs a measured caproate value.")
    return df, genus_cols


def present(cols, df):
    return [c for c in cols if c in df.columns]


def scale(X_train, X_test):
    """Standardize the numeric predictors on the training rows only. Blank cells
    stay blank; categorical columns and the target are never scaled."""
    X_train, X_test = X_train.copy(), X_test.copy()
    num = [c for c in X_train.columns if c not in CAT_COLS]
    if num:
        sc = StandardScaler().fit(X_train[num])
        X_train[num] = sc.transform(X_train[num])
        X_test[num] = sc.transform(X_test[num])
    return X_train, X_test


def r2_rmse(y_true, y_pred):
    return r2_score(y_true, y_pred), root_mean_squared_error(y_true, y_pred)


def training_mean_r2_rmse(y_train, y_test):
    """The simple benchmark (review item 5): predict the training rows' mean caproate."""
    return r2_rmse(y_test, np.full(len(y_test), float(np.mean(y_train))))


def sampling_units(df):
    """One id per sampling unit: the same-day replicate samples of
    REPLICATE_STUDIES share one id; every other row is its own unit. Checks
    that each replicate unit carries a single caproate value."""
    rep = df['Paper_ID'].isin(REPLICATE_STUDIES)
    day_id = df['Paper_ID'] + '|' + df['BIOREACTOR'] + '|' + df['DAY'].astype(str)
    units = pd.Series(np.where(rep, day_id, 'row' + df.index.astype(str)), index=df.index)
    n_values = df.groupby(units)[TARGET].nunique()
    if (n_values > 1).any():
        raise SystemExit(f"Sampling units with more than one caproate value: {n_values[n_values > 1].index.tolist()}")
    return units


def spearman_by_unit(y_true, y_pred, units):
    """Spearman rho between measured and predicted caproate on sampling-unit
    means, so same-day replicates count once. Returns (rho, number of units)."""
    t = pd.DataFrame({'y': np.asarray(y_true, float), 'p': np.asarray(y_pred, float),
                      'u': np.asarray(units)}).groupby('u')[['y', 'p']].mean()
    if len(t) < 3 or t['p'].std() == 0 or t['y'].std() == 0:
        return np.nan, len(t)
    return spearmanr(t['y'], t['p']).statistic, len(t)


def save_table(table, name):
    """Write a results table into the run's dated results folder, when one is set."""
    out_dir = os.environ.get("BIOTWIN_RESULTS_DIR")
    if not out_dir:
        print(f"  (BIOTWIN_RESULTS_DIR not set: {name} not saved)")
        return
    path = os.path.join(out_dir, name)
    if os.path.exists(path):
        raise SystemExit(f"{path} already exists.")
    pd.DataFrame(table).to_csv(path, index=False)
    print(f"  Saved: {path}")
