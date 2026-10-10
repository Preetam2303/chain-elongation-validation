# 09_guild_membership_shap_ranks.py
# -----------------------------------------------------------------------------
# WHICH GENERA, BESIDES THE KNOWN CHAIN ELONGATOR, HOLD A CONSISTENT POSITION?
#
# Caproiciproducens is rank 1 by mean |SHAP| in every fold of every validation
# architecture tested. That is the already-known producer. This script asks the
# question underneath it: which OTHER genera hold a stable rank across folds,
# and do they look like helpers (abundance pushes predicted caproate up) or
# competitors (abundance pushes it down)?
#
# It answers two things at once:
#   (1) the project's original objective -- identify supporting organisms in
#       the guild, not the producer itself;
#   (2) whether "functionally redundant guild" is defensible as a multi-genus
#       claim, or whether the honest scope is strain-level within one genus.
#
# 22 folds across three architectures on the full matrix, matching the paper's
# own schemes, plus one cross-check (fold counts are printed from the run):
#     7  leave-one-study-out            (full matrix; the headline model)
#     5  GroupKFold by physical vessel  (full matrix)
#    10  shuffled 10-fold               (full matrix)
#     1  Illumina-only global model     (Illumina genus matrix, reported separately)
# Every model uses the core predictor set and the shared settings
# (pipeline_settings.py): pH, temperature, HRT, lactate, acetate, ethanol,
# butyrate, propionate, the 2 categoricals and the genera; unmeasured cells
# blank; no subsampling.
#
# For every genus it reports: median rank, how often it lands in the top 5 /
# 10 / 25, the spread of its rank, whether it survives when the Nanopore cohort
# specifically is held out, and the DIRECTION of its effect.
#
# Direction is measured two ways, because they answer different questions:
#   shap_dir  Spearman correlation between a genus's abundance and its own SHAP
#             value, within each fold. Positive = more of this genus pushes the
#             model's caproate prediction UP, holding the rest of the model
#             fixed. This is the model-internal, conditional effect.
#   raw_rho   Spearman correlation between the genus's abundance and measured
#             caproate across all samples, on sampling-unit means (the three
#             same-day replicate samples of Duber 2024 and the Nanopore study
#             count once). Marginal, unconditioned.
#
# Within-study check of the guild claim (review item 7): the Spearman
# correlation of Caproiciproducens abundance with measured caproate inside each
# study separately, on sampling-unit means. A pooled correlation can come from
# differences between studies alone; this shows whether it holds within them.
# A genus with high |SHAP| but shap_dir near zero is being used for splitting
# without a consistent directional effect -- that is NOT a helper claim.
#
# Ties are ranked deterministically (see rank_fold), so genera with EXACTLY zero
# attribution get the same rank on every machine. Genera that are only NEARLY
# tied can still swap places between machines, because SHAP values differ in the
# last decimal places across CPUs; this moves individual counts by about one fold
# (e.g. top-10 in 19 vs 20 of 22) but never changes which genera pass the strict
# or relaxed criteria. In the 'noNano' column, ZERO means the genus received
# exactly zero attribution once the Nanopore cohort was held out -- not a low
# rank, no attribution at all.
# Runtime: about two minutes. Run from inside python/pipeline/.
# -----------------------------------------------------------------------------

import os
import pandas as pd
import numpy as np
import xgboost as xgb
import shap
from scipy.stats import spearmanr, rankdata
from sklearn.model_selection import LeaveOneGroupOut, GroupKFold, KFold
from sklearn.preprocessing import StandardScaler
import warnings

from pipeline_settings import (CAT_COLS, CORE_CHEM, NANOPORE_ID, OPS, TARGET, XGB_PARAMS, load_matrix,
                               present, sampling_units)

warnings.filterwarnings('ignore')

FILE_PATH = os.environ.get("BIOTWIN_MATRIX", "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")
GENUS65_PATH = os.environ.get("BIOTWIN_GENUS", "../../data/historical/BIOTWIN_GENUS_ML_MATRIX.csv")
RESULTS_DIR = os.environ.get("BIOTWIN_RESULTS_DIR", "../../results")
CAT = CAT_COLS
CHEM_OPS = CORE_CHEM + OPS
PRODUCER = 'g__Caproiciproducens'

XP = XGB_PARAMS


def unit_spearman(x, y, units):
    """Spearman rho and n on sampling-unit means."""
    t = pd.DataFrame({'x': np.asarray(x, float), 'y': np.asarray(y, float),
                      'u': np.asarray(units)}).groupby('u')[['x', 'y']].mean()
    if len(t) < 3 or t['x'].std() == 0 or t['y'].std() == 0:
        return np.nan, np.nan, len(t)
    r = spearmanr(t['x'], t['y'])
    return r.statistic, r.pvalue, len(t)


def rank_fold(X_tr, y_tr, params, scale_cols):
    """Fit, compute TreeSHAP, return ({genus: (rank, mean|SHAP|)}, signed-direction dict)."""
    Xs = X_tr.copy()
    sc = StandardScaler()
    Xs[scale_cols] = sc.fit_transform(Xs[scale_cols])
    model = xgb.XGBRegressor(**params).fit(Xs, y_tr)
    sv = shap.TreeExplainer(model).shap_values(Xs)

    cols = list(Xs.columns)
    mean_abs = np.abs(sv).mean(axis=0)
    # Rank genera only. Ties (most often a block of genera with EXACTLY zero
    # attribution) share one averaged rank, so output is identical on every
    # platform. A plain sort breaks ties in platform-dependent order, which made
    # the rank of zero-attribution genera meaningless and non-reproducible.
    gidx = [j for j, c in enumerate(cols) if str(c).startswith('g__')]
    gnames = [cols[j] for j in gidx]
    gvals = mean_abs[gidx]
    granks = rankdata(-gvals, method='average')
    genera = {g: (float(r), float(v)) for g, r, v in zip(gnames, granks, gvals)}

    direction = {}
    for j, c in enumerate(cols):
        if not str(c).startswith('g__'):
            continue
        abund = X_tr[c].values
        shapv = sv[:, j]
        if np.std(abund) == 0 or np.std(shapv) == 0:
            direction[c] = np.nan
        else:
            direction[c] = spearmanr(abund, shapv).statistic
    return genera, direction


def main():
    df, genus_cols = load_matrix(FILE_PATH)
    vco = present(CHEM_OPS, df)
    units = sampling_units(df)

    y = df[TARGET]
    papers = df['Paper_ID']
    df['VID'] = papers + "_" + df['BIOREACTOR']

    records = []   # one row per (fold, genus)
    fold_meta = []

    # ---- 1. LOSO by study (7 folds) ----
    X = df[vco + CAT + genus_cols]
    for tr, te in LeaveOneGroupOut().split(X, y, groups=papers):
        held = papers.iloc[te].unique()[0]
        genera, direction = rank_fold(X.iloc[tr], y.iloc[tr], XP, vco + genus_cols)
        fold_meta.append(('LOSO', f'held out {held}'))
        for g, (r, v) in genera.items():
            records.append({'arch': 'LOSO', 'fold': f'LOSO:{held}', 'genus': g,
                            'rank': r, 'shap': v, 'dir': direction.get(g, np.nan)})

    # ---- 2. GroupKFold by vessel (5 folds) ----
    for i, (tr, te) in enumerate(GroupKFold(n_splits=5).split(X, y, groups=df['VID']), 1):
        genera, direction = rank_fold(X.iloc[tr], y.iloc[tr], XP, vco + genus_cols)
        fold_meta.append(('GroupKFold', f'fold {i}'))
        for g, (r, v) in genera.items():
            records.append({'arch': 'GroupKFold', 'fold': f'GKF:{i}', 'genus': g,
                            'rank': r, 'shap': v, 'dir': direction.get(g, np.nan)})

    # ---- 3. Shuffled 10-fold (10 folds) ----
    for i, (tr, te) in enumerate(KFold(n_splits=10, shuffle=True, random_state=42).split(X), 1):
        genera, direction = rank_fold(X.iloc[tr], y.iloc[tr], XP, vco + genus_cols)
        fold_meta.append(('Shuffled10', f'fold {i}'))
        for g, (r, v) in genera.items():
            records.append({'arch': 'Shuffled10', 'fold': f'S10:{i}', 'genus': g,
                            'rank': r, 'shap': v, 'dir': direction.get(g, np.nan)})

    rec = pd.DataFrame(records)
    n_folds = rec['fold'].nunique()

    # marginal correlation with measured caproate, all samples, on sampling-unit means
    raw_rho = {}
    for g in genus_cols:
        raw_rho[g] = unit_spearman(df[g], y, units)[0]

    # ---- summary ----
    rows = []
    for g, sub in rec.groupby('genus'):
        ranks = sub['rank']
        loso_n = sub[sub['arch'] == 'LOSO'].shape[0]
        nano_row = sub[sub['fold'] == f'LOSO:{NANOPORE_ID}']
        rows.append({
            'genus': g,
            'median_rank': ranks.median(),
            'best': ranks.min(), 'worst': ranks.max(),
            'top5': int((ranks <= 5).sum()), 'top10': int((ranks <= 10).sum()),
            'top25': int((ranks <= 25).sum()),
            'loso_top10': int((sub[sub['arch'] == 'LOSO']['rank'] <= 10).sum()),
            'gkf_top10': int((sub[sub['arch'] == 'GroupKFold']['rank'] <= 10).sum()),
            's10_top10': int((sub[sub['arch'] == 'Shuffled10']['rank'] <= 10).sum()),
            'rank_no_nano': (nano_row['rank'].iloc[0] if len(nano_row) else np.nan),
            'zero_when_no_nano': bool(len(nano_row) and nano_row['shap'].iloc[0] == 0),
            'zero_attr_folds': int((sub['shap'] == 0).sum()),
            'shap_dir': sub['dir'].mean(),
            'raw_rho': raw_rho.get(g, np.nan),
        })
    s = pd.DataFrame(rows).sort_values(['top10', 'median_rank'], ascending=[False, True])

    print("=" * 110)
    print(f"GUILD MEMBERSHIP BY SHAP RANK  --  {n_folds} folds on the {len(df)}-sample, {len(genus_cols)}-genus matrix")
    print("  LOSO = 7 folds | GroupKFold-by-vessel = 5 | Shuffled 10-fold = 10")
    print("=" * 110)
    print(f"{'genus':<34s}{'med':>5s}{'best':>5s}{'wrst':>5s}{'top5':>6s}{'top10':>7s}"
          f"{'top25':>7s}{'LOSO':>6s}{'GKF':>5s}{'S10':>5s}{'noNano':>7s}"
          f"{'shapDir':>9s}{'rawRho':>8s}")
    print("-" * 110)
    for _, r in s.head(20).iterrows():
        nn = ' ZERO' if r['zero_when_no_nano'] else f"{r['rank_no_nano']:5.0f}"
        sd = '    n/a' if np.isnan(r['shap_dir']) else f"{r['shap_dir']:7.3f}"
        rr = '   n/a' if np.isnan(r['raw_rho']) else f"{r['raw_rho']:6.3f}"
        print(f"{r['genus']:<34s}{r['median_rank']:>5.0f}{r['best']:>5.0f}{r['worst']:>5.0f}"
              f"{r['top5']:>4d}/{n_folds:<2d}{r['top10']:>5d}/{n_folds:<2d}"
              f"{r['top25']:>5d}/{n_folds:<2d}{r['loso_top10']:>4d}/7"
              f"{r['gkf_top10']:>3d}/5{r['s10_top10']:>3d}/10{nn:>8s}{sd:>9s}{rr:>8s}")

    # ---- who occupies rank 2 ----
    print("\n" + "=" * 110)
    print("RANK-2 OCCUPANCY  (which genus sits directly below the known producer)")
    print("=" * 110)
    r2 = rec[rec['rank'] == 2].groupby('genus').size().sort_values(ascending=False)
    for g, n in r2.items():
        print(f"  {g:<40s} rank 2 in {n}/{n_folds} folds")

    # ---- strict guild criterion ----
    print("\n" + "=" * 110)
    print("STRICT CRITERION: top-10 in EVERY fold of EVERY architecture")
    print("=" * 110)
    strict = s[(s['loso_top10'] == 7) & (s['gkf_top10'] == 5) & (s['s10_top10'] == 10)]
    if len(strict) == 0:
        print("  NONE. No genus besides the rank-1 producer holds top-10 universally.")
    for _, r in strict.iterrows():
        print(f"  {r['genus']:<40s} median rank {r['median_rank']:.0f}, "
              f"shap_dir {r['shap_dir']:+.3f}, raw_rho {r['raw_rho']:+.3f}")

    print("\n  Relaxed (top-10 in >=80% of folds):")
    rel = s[s['top10'] >= 0.8 * n_folds]
    for _, r in rel.iterrows():
        print(f"  {r['genus']:<40s} top10 {int(r['top10'])}/{n_folds}, "
              f"median {r['median_rank']:.0f}, shap_dir {r['shap_dir']:+.3f}, "
              f"raw_rho {r['raw_rho']:+.3f}")

    # ---- within-study correlation of the producer with caproate ----
    print("\n" + "=" * 110)
    print(f"WITHIN-STUDY CHECK: Spearman rho of {PRODUCER} abundance with measured caproate,")
    print("inside each study, on sampling-unit means (same-day replicates count once)")
    print("=" * 110)
    within = []
    if PRODUCER in df.columns:
        for study, d in df.groupby('Paper_ID'):
            rho, p, n = unit_spearman(d[PRODUCER], d[TARGET], units[d.index])
            within.append({'study': study, 'n_rows': len(d), 'n_sampling_units': n,
                           'mean_rel_abundance_pct': 100 * d[PRODUCER].mean(), 'spearman_rho': rho, 'p_value': p})
            rr = '   n/a' if np.isnan(rho) else f"{rho:+.3f}"
            pp = '' if np.isnan(p) else f" (p = {p:.3f})"
            print(f"  {study:<34s} n = {n:3d} units ({len(d):3d} rows) | mean abundance "
                  f"{100 * d[PRODUCER].mean():5.1f}% | rho {rr}{pp}")
        rho, p, n = unit_spearman(df[PRODUCER], df[TARGET], units)
        print(f"  {'All studies pooled':<34s} n = {n:3d} units ({len(df):3d} rows) | rho {rho:+.3f} (p = {p:.3g})")
        within.append({'study': 'All studies pooled', 'n_rows': len(df), 'n_sampling_units': n,
                       'mean_rel_abundance_pct': 100 * df[PRODUCER].mean(), 'spearman_rho': rho, 'p_value': p})
    else:
        print(f"  {PRODUCER} is not in the matrix.")

    # ---- Illumina-only genus universe, reported separately ----
    try:
        gm, g65 = load_matrix(GENUS65_PATH)
        v65 = present(CHEM_OPS, gm)
        genera65, _ = rank_fold(gm[v65 + CAT + g65], gm[TARGET], XP, v65 + g65)
        print("\n" + "=" * 110)
        print(f"CROSS-CHECK: Illumina-only {len(g65)}-genus universe (single global model)")
        print("=" * 110)
        top65 = sorted(genera65.items(), key=lambda kv: kv[1][0])[:10]
        for g, (r, v) in top65:
            print(f"  {r:4.0f}. {g}")
    except FileNotFoundError:
        print("\n  [Illumina genus matrix not found -- cross-check skipped]")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    for table, name in [(s, "guild_membership_shap_ranks.csv"),
                        (pd.DataFrame(within), "guild_producer_within_study_correlation.csv")]:
        path = os.path.join(RESULTS_DIR, name)
        if os.path.exists(path):
            raise SystemExit(f"{path} already exists.")
        table.to_csv(path, index=False)
        print("\n  Table written to " + path)


if __name__ == "__main__":
    main()
