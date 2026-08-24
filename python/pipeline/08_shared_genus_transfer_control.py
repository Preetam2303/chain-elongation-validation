# 08_shared_genus_transfer_control.py
# -----------------------------------------------------------------------------
# REVIEWER CONTROL for Section 3.6.
#
# The concern: only 19 of the 170 genera in the final matrix are detected on
# BOTH platforms. At test time an Illumina-trained model sees 46 genera forced
# to zero and 105 genera it has no training signal for. So the cross-platform
# failure could be a feature-space artifact rather than an ecological one.
#
# This script tests that directly by rerunning the transfer restricted to the
# shared feature subspace, plus the ablations needed to interpret the answer:
#
#   S0  chemistry + operational only          (no genera at all -- the floor)
#   S1  19 shared genera ONLY                 (no chemistry)
#   S2  all 170 genera ONLY                   (no chemistry)
#   S3  19 shared genera + chemistry + ops    <-- THE KEY TEST
#   S4  18 shared genera (Unclassified dropped) + chemistry + ops
#   S5  all 170 genera + chemistry + ops      (the reported baseline)
#
# Each genus-containing scope is run under two normalizations:
#   [global] relative abundance over all 170 genera, then subset
#   [renorm] relative abundance recomputed WITHIN the retained subset
#
# It also runs two interpretive guards:
#   (a) within-Illumina LOSO at the shared scope -- can 19 genera predict
#       caproate at all, even without crossing platforms?
#   (b) the reverse transfer, Nanopore -> Illumina.
#
# Model configuration is identical to 06_cross_platform_transfer.py, including
# NO subsample/colsample_bytree, so every number here is exactly reproducible
# across platforms (see README's reproducibility note).
#
# Runtime: about a minute. Run from inside python/pipeline/.
# -----------------------------------------------------------------------------

import pandas as pd
import numpy as np
import xgboost as xgb
from scipy.stats import spearmanr
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import root_mean_squared_error, r2_score, mean_absolute_error
import warnings
warnings.filterwarnings('ignore')

FILE_PATH = "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"
NANOPORE_ID = 'Hanna_2025'
TARGET = 'Caproate'
OP_AND_SUBSTRATES = ['PH', 'TEMP', 'HRT', 'Lactate', 'Acetate', 'Ethanol']

XGB_PARAMS = dict(n_estimators=150, learning_rate=0.05, max_depth=3, random_state=42)


def load():
    df = pd.read_csv(FILE_PATH)
    df['Paper_ID'] = df['Paper_ID'].astype(str).str.strip()
    genus_cols = sorted([c for c in df.columns if str(c).startswith('g__')])
    df[genus_cols] = df[genus_cols].fillna(0)
    return df, genus_cols


def evaluate(df, genus_subset, chem_ops, train_mask, renormalize, label):
    """Fit on train_mask rows, test on the rest. Returns a result dict."""
    work = df.copy()
    if genus_subset:
        if renormalize:
            sub_sums = work[genus_subset].sum(axis=1).replace(0, 1)
            work[genus_subset] = work[genus_subset].div(sub_sums, axis=0)
        else:
            all_g = sorted([c for c in work.columns if str(c).startswith('g__')])
            row_sums = work[all_g].sum(axis=1).replace(0, 1)
            work[all_g] = work[all_g].div(row_sums, axis=0)

    features = sorted(list(set(genus_subset + [c for c in chem_ops if c in work.columns])))
    if not features:
        return None

    tr, te = work[train_mask], work[~train_mask]
    X_tr = tr[features].apply(pd.to_numeric, errors='coerce').fillna(0)
    y_tr = pd.to_numeric(tr[TARGET], errors='coerce').fillna(0)
    X_te = te[features].apply(pd.to_numeric, errors='coerce').fillna(0)
    y_te = pd.to_numeric(te[TARGET], errors='coerce').fillna(0)

    model = xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr, y_tr)
    pred = model.predict(X_te)
    rho = np.nan if np.std(pred) == 0 else spearmanr(y_te, pred).statistic

    return {'label': label, 'n_feat': len(features),
            'r2': r2_score(y_te, pred), 'rmse': root_mean_squared_error(y_te, pred),
            'mae': mean_absolute_error(y_te, pred), 'rho': rho}


def loso_within_illumina(df, genus_subset, chem_ops, label):
    """LOSO across the six Illumina studies at a given feature scope."""
    work = df[df['Paper_ID'] != NANOPORE_ID].copy()
    all_g = sorted([c for c in work.columns if str(c).startswith('g__')])
    row_sums = work[all_g].sum(axis=1).replace(0, 1)
    work[all_g] = work[all_g].div(row_sums, axis=0)

    features = sorted(list(set(genus_subset + [c for c in chem_ops if c in work.columns])))
    X = work[features].apply(pd.to_numeric, errors='coerce').fillna(0)
    y = pd.to_numeric(work[TARGET], errors='coerce').fillna(0)
    groups = work['Paper_ID']

    scores = []
    for tr_i, te_i in LeaveOneGroupOut().split(X, y, groups=groups):
        if len(te_i) < 2:
            continue
        # Training-fold-only scaling, matching the convention used throughout
        # this pipeline (Methods 2.6).
        X_tr, X_te = X.iloc[tr_i].copy(), X.iloc[te_i].copy()
        sc = StandardScaler()
        X_tr[features] = sc.fit_transform(X_tr[features])
        X_te[features] = sc.transform(X_te[features])
        m = xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr, y.iloc[tr_i])
        scores.append(r2_score(y.iloc[te_i], m.predict(X_te)))
    return {'label': label, 'n_feat': len(features),
            'r2': float(np.mean(scores)), 'sd': float(np.std(scores)),
            'folds': [round(s, 3) for s in scores]}


def row(r):
    rho = "   n/a" if (r['rho'] is None or np.isnan(r['rho'])) else f"{r['rho']:6.3f}"
    return (f"  {r['label']:<46s} {r['n_feat']:>4d}  {r['r2']:>8.3f}  "
            f"{r['rmse']:>7.1f}  {r['mae']:>7.1f}  {rho}")


def main():
    df, genus_cols = load()
    illumina = df['Paper_ID'] != NANOPORE_ID

    ill_present = set(c for c in genus_cols if (df.loc[illumina, c] > 0).any())
    nan_present = set(c for c in genus_cols if (df.loc[~illumina, c] > 0).any())
    shared = sorted(ill_present & nan_present)
    shared_no_unc = [c for c in shared if c != 'g__Unclassified']

    print("=" * 96)
    print("SHARED FEATURE SPACE DIAGNOSTIC")
    print("=" * 96)
    print(f"  Total genus columns in matrix        : {len(genus_cols)}")
    print(f"  Detected on Illumina                 : {len(ill_present)}")
    print(f"  Detected on Nanopore                 : {len(nan_present)}")
    print(f"  SHARED (detected on both)            : {len(shared)}")
    print(f"  Illumina-only (forced to 0 at test)  : {len(ill_present - nan_present)}")
    print(f"  Nanopore-only (no training signal)   : {len(nan_present - ill_present)}")

    rel = df[genus_cols].div(df[genus_cols].sum(axis=1).replace(0, 1), axis=0)
    print(f"\n  Compositional mass carried by the shared set:")
    print(f"    Illumina samples : {100 * rel.loc[illumina, shared].sum(axis=1).mean():5.1f}%")
    print(f"    Nanopore samples : {100 * rel.loc[~illumina, shared].sum(axis=1).mean():5.1f}%")
    if 'g__Unclassified' in shared:
        print(f"  Of which g__Unclassified alone:")
        print(f"    Illumina : {100 * rel.loc[illumina, 'g__Unclassified'].mean():5.1f}%"
              f"    Nanopore : {100 * rel.loc[~illumina, 'g__Unclassified'].mean():5.1f}%")

    print(f"\n  The {len(shared)} shared genera:")
    for g in shared:
        print(f"    {g}")

    header = (f"\n  {'Feature scope':<46s} {'nFeat':>4s}  {'R2':>8s}  "
              f"{'RMSE':>7s}  {'MAE':>7s}  {'rho':>6s}")

    for renorm, tag in [(False, "global rel. abundance over all 170, then subset"),
                        (True,  "rel. abundance RENORMALIZED within retained subset")]:
        print("\n" + "=" * 96)
        print(f"ILLUMINA -> NANOPORE TRANSFER   [{tag}]")
        print("=" * 96 + header)
        # 2 x 3 design: {chem only, shared genera, all genera} x {with, without chem}
        scopes = [
            ([],            [],                 "S0  chemistry + ops only"),
            (shared,        [],                 f"S1  {len(shared)} shared genera ONLY (no chemistry)"),
            (genus_cols,    [],                 "S2  all 170 genera ONLY (no chemistry)"),
            (shared,        OP_AND_SUBSTRATES,  f"S3  {len(shared)} shared genera + chemistry + ops   <-- KEY"),
            (shared_no_unc, OP_AND_SUBSTRATES,  f"S4  {len(shared_no_unc)} shared genera (no Unclassified) + chem + ops"),
            (genus_cols,    OP_AND_SUBSTRATES,  "S5  all 170 genera + chemistry + ops  [BASELINE]"),
        ]
        for gs, chem, label in scopes:
            if not gs and renorm:
                continue  # S0 has no genera; identical under both normalizations
            r = evaluate(df, gs, chem, illumina, renorm, label)
            if r:
                print(row(r))

    print("\n" + "=" * 96)
    print("GUARD (a): WITHIN-ILLUMINA LOSO -- can this feature scope predict at all?")
    print("=" * 96)
    print(f"  {'Feature scope':<46s} {'nFeat':>4s}  {'R2 mean':>8s}  {'SD':>7s}   folds")
    guard_arms = [
        (shared,              OP_AND_SUBSTRATES, f"{len(shared)} shared genera + chemistry + ops"),
        (sorted(ill_present), OP_AND_SUBSTRATES, "65 Illumina-detected genera + chem + ops"),
        (shared,              [],                f"{len(shared)} shared genera ONLY"),
        (sorted(ill_present), [],                "65 Illumina-detected genera ONLY"),
    ]
    for gs, chem, label in guard_arms:
        g = loso_within_illumina(df, gs, chem, label)
        print(f"  {g['label']:<46s} {g['n_feat']:>4d}  {g['r2']:>8.3f}  {g['sd']:>7.3f}   {g['folds']}")

    print("\n" + "=" * 96)
    print("GUARD (b): REVERSE TRANSFER, NANOPORE -> ILLUMINA")
    print("=" * 96 + header)
    for gs, label, chem in [(shared, f"{len(shared)} shared genera + chemistry + ops", OP_AND_SUBSTRATES),
                            (genus_cols, "all 170 genera + chemistry + ops", OP_AND_SUBSTRATES)]:
        r = evaluate(df, gs, chem, ~illumina, False, label)
        if r:
            print(row(r))

    print("\n" + "=" * 96)
    print("SANITY GATE")
    print("=" * 96)
    base = evaluate(df, genus_cols, OP_AND_SUBSTRATES, illumina, False, "baseline")
    ok = abs(base['r2'] - (-0.170)) < 0.002 and abs(base['rmse'] - 145.6) < 0.2
    print(f"  S5 baseline reproduces Table 3b (-0.170 / 145.6 mM C)? "
          f"got {base['r2']:.3f} / {base['rmse']:.1f}  ->  {'PASS' if ok else 'FAIL'}")
    if not ok:
        print("  If this FAILS, something upstream differs and the numbers above")
        print("  should not be interpreted until it is resolved.")


if __name__ == "__main__":
    main()
