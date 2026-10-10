# 08_shared_genus_transfer_control.py
# -----------------------------------------------------------------------------
# REVIEWER CONTROL for Section 3.6 (Table S6).
#
# The concern: only a minority of the genera in the final matrix are detected
# on BOTH platforms. At test time an Illumina-trained model sees the
# Illumina-only genera forced to zero, and the Nanopore-only genera carry no
# training signal. So the cross-platform failure could be a feature-space
# artifact rather than an ecological one. (The counts are computed from the
# matrix and printed below; the published matrix had 19 shared of 170.)
#
# This script tests that directly by rerunning the transfer restricted to the
# shared feature subspace, plus the ablations needed to interpret the answer:
#
#   S0  core chemistry + operational + categoricals only (no genera -- the floor)
#   S1  shared genera ONLY                    (no chemistry)
#   S2  all genera ONLY                       (no chemistry)
#   S3  shared genera + core set              <-- THE KEY TEST
#   S4  shared genera (Unclassified dropped) + core set
#   S5  all genera + core set                 (= 06's relative-abundance transfer)
#
# Each genus-containing scope is run under two normalizations:
#   [global] relative abundance over all genera, then subset
#   [renorm] relative abundance recomputed WITHIN the retained subset
#
# It also runs two interpretive guards:
#   (a) within-Illumina LOSO at the shared scope -- can the shared genera
#       predict caproate at all, even without crossing platforms?
#   (b) the reverse transfer, Nanopore -> Illumina.
#
# Core set and model settings come from pipeline_settings.py (no subsampling),
# so every number here is exactly reproducible across platforms. Spearman rho
# is computed on sampling-unit means (same-day replicates count once).
#
# Runtime: about a minute. Run from inside python/pipeline/.
# -----------------------------------------------------------------------------

import os
import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import mean_absolute_error
import warnings

from pipeline_settings import (CAT_COLS, CORE_CHEM, NANOPORE_ID, OPS, TARGET, XGB_PARAMS, load_matrix,
                               r2_rmse, sampling_units, scale, spearman_by_unit)

warnings.filterwarnings('ignore')

FILE_PATH = os.environ.get("BIOTWIN_MATRIX", "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")
OP_AND_SUBSTRATES = OPS + CORE_CHEM + CAT_COLS


def load():
    df, genus_cols = load_matrix(FILE_PATH, relative=False)  # normalized per scope below
    return df, genus_cols, sampling_units(df)


def evaluate(df, units, genus_subset, chem_ops, train_mask, renormalize, label):
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
    model = xgb.XGBRegressor(**XGB_PARAMS).fit(tr[features], tr[TARGET])
    pred = model.predict(te[features])
    r2, rmse = r2_rmse(te[TARGET], pred)
    rho, _ = spearman_by_unit(te[TARGET], pred, units[te.index])

    return {'label': label, 'n_feat': len(features), 'r2': r2, 'rmse': rmse,
            'mae': mean_absolute_error(te[TARGET], pred), 'rho': rho}


def loso_within_illumina(df, genus_subset, chem_ops, label):
    """LOSO across the six Illumina studies at a given feature scope."""
    work = df[df['Paper_ID'] != NANOPORE_ID].copy()
    all_g = sorted([c for c in work.columns if str(c).startswith('g__')])
    row_sums = work[all_g].sum(axis=1).replace(0, 1)
    work[all_g] = work[all_g].div(row_sums, axis=0)

    features = sorted(list(set(genus_subset + [c for c in chem_ops if c in work.columns])))
    X, y, groups = work[features], work[TARGET], work['Paper_ID']

    scores = []
    for tr_i, te_i in LeaveOneGroupOut().split(X, y, groups=groups):
        if len(te_i) < 2:
            continue
        # Training-fold-only scaling, matching the convention used throughout
        # this pipeline (Methods 2.6).
        X_tr, X_te = scale(X.iloc[tr_i], X.iloc[te_i])
        m = xgb.XGBRegressor(**XGB_PARAMS).fit(X_tr, y.iloc[tr_i])
        scores.append(r2_rmse(y.iloc[te_i], m.predict(X_te))[0])
    return {'label': label, 'n_feat': len(features),
            'r2': float(np.mean(scores)), 'sd': float(np.std(scores)),
            'folds': [round(s, 3) for s in scores]}


def row(r):
    rho = "   n/a" if (r['rho'] is None or np.isnan(r['rho'])) else f"{r['rho']:6.3f}"
    return (f"  {r['label']:<46s} {r['n_feat']:>4d}  {r['r2']:>8.3f}  "
            f"{r['rmse']:>7.1f}  {r['mae']:>7.1f}  {rho}")


def main():
    df, genus_cols, units = load()
    illumina = df['Paper_ID'] != NANOPORE_ID
    n_all = len(genus_cols)

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

    for renorm, tag in [(False, f"global rel. abundance over all {n_all}, then subset"),
                        (True,  "rel. abundance RENORMALIZED within retained subset")]:
        print("\n" + "=" * 96)
        print(f"ILLUMINA -> NANOPORE TRANSFER   [{tag}]")
        print("=" * 96 + header)
        # 2 x 3 design: {chem only, shared genera, all genera} x {with, without chem}
        scopes = [
            ([],            OP_AND_SUBSTRATES,  "S0  core chemistry + ops + categoricals only"),
            (shared,        [],                 f"S1  {len(shared)} shared genera ONLY (no chemistry)"),
            (genus_cols,    [],                 f"S2  all {n_all} genera ONLY (no chemistry)"),
            (shared,        OP_AND_SUBSTRATES,  f"S3  {len(shared)} shared genera + chemistry + ops   <-- KEY"),
            (shared_no_unc, OP_AND_SUBSTRATES,  f"S4  {len(shared_no_unc)} shared genera (no Unclassified) + chem + ops"),
            (genus_cols,    OP_AND_SUBSTRATES,  f"S5  all {n_all} genera + chemistry + ops  [BASELINE]"),
        ]
        for gs, chem, label in scopes:
            if not gs and renorm:
                continue  # S0 has no genera; identical under both normalizations
            r = evaluate(df, units, gs, chem, illumina, renorm, label)
            if r:
                print(row(r))

    print("\n" + "=" * 96)
    print("GUARD (a): WITHIN-ILLUMINA LOSO -- can this feature scope predict at all?")
    print("=" * 96)
    print(f"  {'Feature scope':<46s} {'nFeat':>4s}  {'R2 mean':>8s}  {'SD':>7s}   folds")
    guard_arms = [
        (shared,              OP_AND_SUBSTRATES, f"{len(shared)} shared genera + chemistry + ops"),
        (sorted(ill_present), OP_AND_SUBSTRATES, f"{len(ill_present)} Illumina-detected genera + chem + ops"),
        (shared,              [],                f"{len(shared)} shared genera ONLY"),
        (sorted(ill_present), [],                f"{len(ill_present)} Illumina-detected genera ONLY"),
    ]
    for gs, chem, label in guard_arms:
        g = loso_within_illumina(df, gs, chem, label)
        print(f"  {g['label']:<46s} {g['n_feat']:>4d}  {g['r2']:>8.3f}  {g['sd']:>7.3f}   {g['folds']}")

    print("\n" + "=" * 96)
    print("GUARD (b): REVERSE TRANSFER, NANOPORE -> ILLUMINA")
    print("=" * 96 + header)
    for gs, label, chem in [(shared, f"{len(shared)} shared genera + chemistry + ops", OP_AND_SUBSTRATES),
                            (genus_cols, f"all {n_all} genera + chemistry + ops", OP_AND_SUBSTRATES)]:
        r = evaluate(df, units, gs, chem, ~illumina, False, label)
        if r:
            print(row(r))

    print("\n" + "=" * 96)
    print("CONSISTENCY CHECK")
    print("=" * 96)
    base = evaluate(df, units, genus_cols, OP_AND_SUBSTRATES, illumina, False, "baseline")
    print(f"  S5 (all genera + core set, relative abundance) = {base['r2']:.3f} / {base['rmse']:.1f} mM C.")
    print("  This must equal 06_cross_platform_transfer.py's relative-abundance transfer;")
    print("  if it does not, the two scripts are not using the same data or settings.")
    ref = os.path.join(os.environ.get("BIOTWIN_RESULTS_DIR", ""), "06_platform_transfer.csv")
    if os.path.exists(ref):
        t06 = pd.read_csv(ref)
        r06 = t06.loc[t06['scheme'] == 'Platform transfer, relative abundance', 'r2'].iloc[0]
        ok = abs(base['r2'] - r06) < 1e-9
        print(f"  06 gave {r06:.3f}  ->  {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
