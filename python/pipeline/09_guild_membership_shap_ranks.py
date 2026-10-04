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
# 23 folds across four architectures, matching the paper's own schemes:
#     7  leave-one-study-out            (170-genus matrix)
#     5  GroupKFold by physical vessel  (170-genus matrix)
#    10  shuffled 10-fold               (170-genus matrix)
#     1  Illumina-only global model     (65-genus matrix, reported separately)
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
#             caproate across all 184 samples. Marginal, unconditioned.
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

import pandas as pd
import numpy as np
import xgboost as xgb
import shap
from scipy.stats import spearmanr, rankdata
from sklearn.model_selection import LeaveOneGroupOut, GroupKFold, KFold
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

FILE_PATH = "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"
GENUS65_PATH = "../../data/historical/BIOTWIN_GENUS_ML_MATRIX.csv"
NANOPORE_ID = 'Hanna_2025'
TARGET = 'Caproate'
CAT = ['Feed_Complexity', 'Primary_Carbon_Signature']
CHEM_OPS = ['Lactate', 'Acetate', 'Butyrate', 'Ethanol', 'Propionate', 'PH', 'TEMP', 'HRT']

XP = dict(n_estimators=150, learning_rate=0.05, max_depth=3, random_state=42,
          n_jobs=-1, enable_categorical=True, tree_method='hist')
XP_NOCAT = {k: v for k, v in XP.items() if k != 'enable_categorical'}


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
    df = pd.read_csv(FILE_PATH)
    df['Paper_ID'] = df['Paper_ID'].astype(str).str.strip()
    genus_cols = sorted([c for c in df.columns if str(c).startswith('g__')])
    df[genus_cols] = df[genus_cols].fillna(0)
    df[genus_cols] = df[genus_cols].div(df[genus_cols].sum(axis=1).replace(0, 1), axis=0)
    for c in CAT:
        df[c] = df[c].astype('category')
    vco = [c for c in CHEM_OPS if c in df.columns]
    for c in vco + genus_cols + [TARGET]:
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0)

    y = df[TARGET]
    papers = df['Paper_ID']
    df['VID'] = papers + "_" + df['BIOREACTOR'].astype(str).str.strip()

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
    Xv = df[vco + genus_cols]
    for i, (tr, te) in enumerate(GroupKFold(n_splits=5).split(Xv, y, groups=df['VID']), 1):
        genera, direction = rank_fold(Xv.iloc[tr], y.iloc[tr], XP_NOCAT, vco + genus_cols)
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

    # marginal correlation with measured caproate, full 184 samples
    raw_rho = {}
    for g in genus_cols:
        raw_rho[g] = (np.nan if df[g].std() == 0
                      else spearmanr(df[g], y).statistic)

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
    print(f"GUILD MEMBERSHIP BY SHAP RANK  --  {n_folds} folds on the 184-sample, 170-genus matrix")
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

    # ---- Illumina-only 65-genus universe, reported separately ----
    try:
        gm = pd.read_csv(GENUS65_PATH)
        g65 = sorted([c for c in gm.columns if str(c).startswith('g__')])
        v65 = [c for c in CHEM_OPS if c in gm.columns]
        Xg = gm[v65 + g65].apply(pd.to_numeric, errors='coerce').fillna(0)
        yg = pd.to_numeric(gm[TARGET], errors='coerce').fillna(0)
        genera65, _ = rank_fold(Xg, yg, XP_NOCAT, v65 + g65)
        print("\n" + "=" * 110)
        print("CROSS-CHECK: Illumina-only 65-genus universe (single global model)")
        print("=" * 110)
        top65 = sorted(genera65.items(), key=lambda kv: kv[1][0])[:10]
        for g, (r, v) in top65:
            print(f"  {r:4.0f}. {g}")
    except FileNotFoundError:
        print("\n  [65-genus matrix not found at data/historical/ -- cross-check skipped]")

    s.to_csv("../../results/guild_membership_shap_ranks.csv", index=False)
    print("\n  Full table written to results/guild_membership_shap_ranks.csv")


if __name__ == "__main__":
    main()
