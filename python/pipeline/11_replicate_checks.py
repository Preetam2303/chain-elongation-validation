# 11_replicate_checks.py
# The two replicate checks planned on 2026-10-07 (project file
# rerun/2026-10-07_replicates/replicates_plan.md, section 5).
#
# Duber 2024 and the Nanopore study sequenced three samples from the same
# reactor on each sampling day, and all three carry that day's single
# chemistry measurement, caproate included (sampling_units() in
# pipeline_settings.py). The main analysis keeps every sample as a row, and
# every grouped scheme keeps a day's three samples together. These checks show
# what keeping them does:
#
#   1. Replicate similarity. Bray-Curtis dissimilarity between the genus
#      profiles (relative abundance) of two samples from the same reactor:
#      from the same day, from consecutive sampling days, and from any two
#      different days. Also: for how many samples the most similar other
#      sample of the study is one of its own same-day replicates.
#   2. One row per sampling day. Each day's three samples are averaged into one
#      row (the genera or ASVs; caproate and the other predictors are the same
#      within a day), and the headline LOSO (02) and the four naive schemes of
#      Table 3a (10) are rerun on those matrices. The predictors, settings and
#      seeds are those of 02 and 10, so the every-sample rows reproduce their
#      numbers.
#   3. The naive schemes with other seeds. Which rows land in the test set
#      changes a naive score a lot, so one seed can make the two versions look
#      closer or further apart than they are. The four naive schemes of 2 are
#      also run with seeds 0-19, at both levels; LOSO has no random split.
#
# Output (in BIOTWIN_RESULTS_DIR): 11_replicate_bray_curtis.csv,
# 11_one_row_per_day.csv and 11_naive_splits_over_seeds.csv.
# Runtime: section 3 fits 20 x 2 x 4 naive schemes and takes the longest.

import os
import warnings
from itertools import combinations

import numpy as np
import pandas as pd
import shap
import xgboost as xgb
from scipy.spatial.distance import braycurtis, pdist, squareform
from sklearn.model_selection import KFold, LeaveOneGroupOut, train_test_split

from pipeline_settings import (CAT_COLS, CORE_CHEM, OPS, REPLICATE_STUDIES, TARGET, XGB_PARAMS, drop_unmodelled,
                               load_matrix, present, r2_rmse, sampling_units, save_table, scale,
                               training_mean_r2_rmse)

warnings.filterwarnings('ignore')

MATRIX_PATH = os.environ.get("BIOTWIN_MATRIX", "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")
PRUNED_PATH = os.environ.get("BIOTWIN_PRUNED", "../../data/historical/BIOTWIN_PRUNED_ML_MATRIX.csv")
CORE_NON_GENUS = OPS + CORE_CHEM + CAT_COLS
LEVELS = ['every sample', 'one row per sampling day']


def load_asv_matrix(path):
    """The ASV-level matrix, loaded as in 10_table3a_naive_splits.py."""
    df = pd.read_csv(path, low_memory=False)
    for c in ['Paper_ID', 'Sample_ID', 'BIOREACTOR']:
        df[c] = df[c].astype(str).str.strip()
    df = drop_unmodelled(df)
    for c in OPS + CORE_CHEM + [TARGET]:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    for c in CAT_COLS:
        df[c] = df[c].astype('category')
    if df[TARGET].isna().any():
        raise SystemExit("Caproate is missing in the ASV-level matrix; every modelled row needs it.")
    return df, [c for c in df.columns if c.startswith('ASV_')]


def one_row_per_day(df, averaged):
    """One row per sampling unit. The averaged columns (genera or ASVs) are
    averaged over the unit's samples; every other column except Sample_ID must
    be the same within a unit (checked), and the unit's Sample_IDs are joined.
    The categorical columns keep their categories, so XGBoost sees the same codes."""
    groups = df.groupby(sampling_units(df).values, sort=False)
    others = [c for c in df.columns if c not in averaged and c != 'Sample_ID']
    varying = [c for c in others if (groups[c].nunique(dropna=False) > 1).any()]
    if varying:
        raise SystemExit(f"Columns that differ between the samples of one sampling day: {varying}")
    out = groups[others].first()
    out[averaged] = groups[averaged].mean()
    out['Sample_ID'] = groups['Sample_ID'].agg('+'.join)
    out = out.reset_index(drop=True)[list(df.columns)]
    for c in CAT_COLS:
        if c in out.columns:
            out[c] = out[c].astype(df[c].dtype)
    return out


# ==============================================================================
# 1. REPLICATE SIMILARITY (Bray-Curtis on relative genus abundance)
# ==============================================================================
def bray_curtis_pairs(df, genus_cols):
    """Every pair of samples from the same reactor of a replicate study, with
    how many sampling days apart they are (0 = the same day). The Nanopore DAY
    labels carry a 'd' prefix ('d15'), so the day number is read from the label."""
    pairs = []
    rep = df[df['Paper_ID'].isin(REPLICATE_STUDIES)]
    for (study, reactor), d in rep.groupby(['Paper_ID', 'BIOREACTOR']):
        day = d['DAY'].astype(str).str.extract(r'(\d+\.?\d*)', expand=False).astype(float).to_numpy()
        if np.isnan(day).any():
            raise SystemExit(f"DAY labels without a day number in {study} {reactor}: {sorted(set(d['DAY']))}")
        day_number = {v: i for i, v in enumerate(sorted(set(day)))}
        X = d[genus_cols].to_numpy(float)
        for i, j in combinations(range(len(d)), 2):
            pairs.append({'study': study, 'reactor': reactor,
                          'days_apart': abs(day_number[day[i]] - day_number[day[j]]),
                          'bray_curtis': braycurtis(X[i], X[j])})
    return pd.DataFrame(pairs)


def nearest_is_replicate(df, genus_cols, study):
    """For each sample of the study: is the most similar other sample of the
    study (lowest Bray-Curtis) one of its own same-day replicates?"""
    d = df[df['Paper_ID'] == study]
    dist = squareform(pdist(d[genus_cols].to_numpy(float), metric='braycurtis'))
    np.fill_diagonal(dist, np.inf)
    unit = (d['BIOREACTOR'] + '|' + d['DAY'].astype(str)).to_numpy()
    nearest = dist.argmin(axis=1)
    return int((unit[nearest] == unit).sum()), len(d)


full, genus_cols = load_matrix(MATRIX_PATH)
pairs = bray_curtis_pairs(full, genus_cols)
kinds = [('same day', pairs['days_apart'] == 0),
         ('next sampling day', pairs['days_apart'] == 1),
         ('any two different days', pairs['days_apart'] >= 1)]

print("--- 1. REPLICATE SIMILARITY: Bray-Curtis dissimilarity of genus profiles (0 = identical, 1 = nothing shared) ---")
print("Pairs of samples from the same reactor")
print("-" * 110)
bc_table = []
for study in REPLICATE_STUDIES:
    in_study = pairs['study'] == study
    k, n = nearest_is_replicate(full, genus_cols, study)
    print(f"{study} ({in_study.sum()} pairs in {pairs.loc[in_study, 'reactor'].nunique()} reactors)")
    for label, sel in kinds:
        v = pairs.loc[in_study & sel, 'bray_curtis'].to_numpy()
        q = np.percentile(v, [0, 25, 50, 75, 100])
        print(f"  {label:24s} {len(v):4d} pairs | median {q[2]:.3f} (middle half {q[1]:.3f}-{q[3]:.3f}) | "
              f"range {q[0]:.3f}-{q[4]:.3f}")
        bc_table.append({'study': study, 'pairs': label, 'n_pairs': len(v), 'median': q[2], 'q25': q[1],
                         'q75': q[3], 'min': q[0], 'max': q[4],
                         'samples_nearest_to_own_replicate': k, 'samples': n})
    print(f"  For {k} of {n} samples, the most similar other sample of the study is a same-day replicate.")
print("-" * 110)
save_table(bc_table, "11_replicate_bray_curtis.csv")


# ==============================================================================
# 2. ONE ROW PER SAMPLING DAY: the headline LOSO and Table 3a's naive schemes
# ==============================================================================
def split_sets(units, fold_of_row):
    """Replicate sets (units with more than one row) whose rows fall in more than one fold."""
    t = pd.DataFrame({'u': units.values, 'f': np.asarray(fold_of_row)})
    multi = t.groupby('u').filter(lambda g: len(g) > 1)
    if multi.empty:
        return 0, 0
    n_folds = multi.groupby('u')['f'].nunique()
    return int((n_folds > 1).sum()), int(len(n_folds))


def loso(data, features):
    """The headline: LOSO with the core set (02_loso_by_study.py, tier 3), for
    XGBoost and the training-mean benchmark."""
    X, y, groups = data[features], data[TARGET], data['Paper_ID']
    res = {'XGBoost': ([], []), 'Training-mean benchmark': ([], [])}
    held_out = []
    fold_of_row = np.zeros(len(data), int)
    for f, (tr, te) in enumerate(LeaveOneGroupOut().split(X, y, groups=groups), 1):
        fold_of_row[te] = f
        if len(te) < 2:
            continue
        X_tr, X_te = scale(X.iloc[tr], X.iloc[te])
        model = xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr, y.iloc[tr])
        for name, (r2, rmse) in [('XGBoost', r2_rmse(y.iloc[te], model.predict(X_te))),
                                 ('Training-mean benchmark', training_mean_r2_rmse(y.iloc[tr], y.iloc[te]))]:
            res[name][0].append(r2); res[name][1].append(rmse)
        held_out.append(groups.iloc[te].iloc[0])
    split, total = split_sets(sampling_units(data), fold_of_row)
    return [{'model': name, 'r2_mean': np.mean(r2s), 'r2_sd': np.std(r2s), 'rmse_mean': np.mean(rmses),
             'rmse_sd': np.std(rmses), 'n': len(data), 'n_features': 0 if name.startswith('Training') else len(features),
             'replicate_sets_split': split, 'replicate_sets': total,
             'folds': '; '.join(f"{h} {r:.3f}" for h, r in zip(held_out, r2s))}
            for name, (r2s, rmses) in res.items()]


def single_split(df, features, seed=42):
    """10_table3a_naive_splits.py's single 80/20 split."""
    units = sampling_units(df)
    idx_train, idx_test = train_test_split(np.arange(len(df)), test_size=0.2, random_state=seed)
    X_tr, X_te = scale(df[features].iloc[idx_train], df[features].iloc[idx_test])
    model = xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr, df[TARGET].iloc[idx_train])
    r2, rmse = r2_rmse(df[TARGET].iloc[idx_test], model.predict(X_te))
    side = np.zeros(len(df), int); side[idx_test] = 1
    split, total = split_sets(units, side)
    return [{'model': 'XGBoost', 'r2_mean': r2, 'r2_sd': np.nan, 'rmse_mean': rmse, 'rmse_sd': np.nan,
             'n': len(df), 'n_features': len(features), 'replicate_sets_split': split, 'replicate_sets': total,
             'folds': ''}]


def kfold(df, features, k, top25=False, seed=42):
    """10_table3a_naive_splits.py's shuffled k-fold CV, optionally with per-fold SHAP top-25 genera."""
    units = sampling_units(df)
    r2s, rmses = [], []
    fold_of_row = np.zeros(len(df), int)
    for f, (tr, te) in enumerate(KFold(k, shuffle=True, random_state=seed).split(df), 1):
        fold_of_row[te] = f
        X_tr, X_te = scale(df[features].iloc[tr], df[features].iloc[te])
        cols = features
        if top25:
            explore = xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr, df[TARGET].iloc[tr])
            sv = np.abs(shap.TreeExplainer(explore).shap_values(X_tr)).mean(axis=0)
            ranked = pd.Series(sv, index=features).sort_values(ascending=False, kind='mergesort')
            top = [c for c in ranked.index if c.startswith('g__')][:25]
            cols = sorted(set([c for c in features if not c.startswith('g__')] + top))
        model = xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr[cols], df[TARGET].iloc[tr])
        r2, rmse = r2_rmse(df[TARGET].iloc[te], model.predict(X_te[cols]))
        r2s.append(r2); rmses.append(rmse)
    split, total = split_sets(units, fold_of_row)
    return [{'model': 'XGBoost', 'r2_mean': np.mean(r2s), 'r2_sd': np.std(r2s), 'rmse_mean': np.mean(rmses),
             'rmse_sd': np.std(rmses), 'n': len(df), 'n_features': len(features),
             'replicate_sets_split': split, 'replicate_sets': total,
             'folds': '; '.join(f"{r:.3f}" for r in r2s)}]


def naive_schemes(asv_m, full_m, seed=42):
    """The four naive schemes of Table 3a (10_table3a_naive_splits.py) on one pair of matrices."""
    d24 = full_m[full_m['Paper_ID'] == 'Duber_2024'].reset_index(drop=True)
    return [
        ("Single random 80/20 split (ASV level)", single_split(asv_m, asv_features, seed)),
        ("Random 5-fold CV (ASV level)", kfold(asv_m, asv_features, 5, seed=seed)),
        ("Shuffled 10-fold CV, per-fold SHAP top-25 genera", kfold(full_m, full_features, 10, top25=True, seed=seed)),
        ("Single-study internal 80/20 split (Duber 2024)", single_split(d24, full_features, seed)),
    ]


asv, asv_cols = load_asv_matrix(PRUNED_PATH)
asv_features = sorted(set(present(CORE_NON_GENUS, asv))) + asv_cols
full_features = sorted(set(present(CORE_NON_GENUS, full) + genus_cols))
matrices = {'every sample': (asv, full),
            'one row per sampling day': (one_row_per_day(asv, asv_cols), one_row_per_day(full, genus_cols))}

table = []
for level in LEVELS:
    asv_m, full_m = matrices[level]
    runs = [("LOSO, core set (the headline)", loso(full_m, full_features))] + naive_schemes(asv_m, full_m)
    for scheme, rows in runs:
        for r in rows:
            table.append({'scheme': scheme, 'rows': level, **r})

print("\n--- 2. ONE ROW PER SAMPLING DAY: each day's replicate samples averaged into one row ---")
print(f"Full matrix: {len(full)} samples -> {len(matrices[LEVELS[1]][1])} rows; "
      f"ASV-level matrix: {len(asv)} samples -> {len(matrices[LEVELS[1]][0])} rows")
print("-" * 110)
print(f"{'':58s}{'every sample':>27s}{'one row per sampling day':>29s}")
res = pd.DataFrame(table)
for (scheme, model), g in res.groupby(['scheme', 'model'], sort=False):
    cells = []
    for level in LEVELS:
        r = g[g['rows'] == level].iloc[0]
        sd = '' if np.isnan(r['r2_sd']) else f" ± {r['r2_sd']:.3f}"
        cells.append(f"R2 {r['r2_mean']:6.3f}{sd:9s} n {r['n']:3d}")
    label = scheme if model == 'XGBoost' else f"{scheme}: {model}"
    print(f"  {label:56s}{cells[0]:>27s}{cells[1]:>29s}")
print("RMSE, mM C:")
for (scheme, model), g in res.groupby(['scheme', 'model'], sort=False):
    vals = [g[g['rows'] == level].iloc[0]['rmse_mean'] for level in LEVELS]
    label = scheme if model == 'XGBoost' else f"{scheme}: {model}"
    print(f"  {label:56s}{vals[0]:>27.1f}{vals[1]:>29.1f}")
headline = res[(res['scheme'].str.contains('headline')) & (res['model'] == 'XGBoost')]
for _, r in headline.iterrows():
    print(f"\nHeadline LOSO, {r['rows']}, R2 per held-out study: {r['folds']}")
print("-" * 110)
save_table(table, "11_one_row_per_day.csv")


# ==============================================================================
# 3. THE NAIVE SCHEMES WITH SEEDS 0-19
# ==============================================================================
SEEDS = range(20)
over_seeds = []
for seed in SEEDS:
    for level in LEVELS:
        for scheme, rows in naive_schemes(*matrices[level], seed=seed):
            over_seeds.append({'seed': seed, 'scheme': scheme, 'rows': level,
                               'r2': rows[0]['r2_mean'], 'rmse': rows[0]['rmse_mean']})

print(f"\n--- 3. THE NAIVE SCHEMES WITH SEEDS {min(SEEDS)}-{max(SEEDS)} (section 2 used seed 42): "
      "R2 median [lowest, highest] ---")
print("-" * 110)
print(f"{'':52s}{'every sample':>20s}{'one row per sampling day':>28s}  lower with one row per day")
seeds_df = pd.DataFrame(over_seeds)
for scheme, g in seeds_df.groupby('scheme', sort=False):
    w = g.pivot(index='seed', columns='rows', values='r2')
    cells = [f"{w[level].median():6.3f} [{w[level].min():.2f}, {w[level].max():.2f}]" for level in LEVELS]
    lower = int((w[LEVELS[1]] < w[LEVELS[0]]).sum())
    print(f"  {scheme:50s}{cells[0]:>20s}{cells[1]:>28s}  in {lower} of {len(w)} seeds")
print("-" * 110)
save_table(over_seeds, "11_naive_splits_over_seeds.csv")
