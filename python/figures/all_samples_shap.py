# all_samples_shap.py
# The one model behind Figures 5 and S5: one XGBoost fit on every modelled
# sample of the corrected matrix (BIOTWIN_MATRIX, from 00b: the Illumina and
# Nanopore samples), with the headline's predictors (the core set: pH,
# temperature, HRT, lactate, acetate, ethanol, butyrate, propionate, the 2
# categorical descriptors, and the genera as relative abundance) and the shared
# settings in pipeline_settings.py (no subsampling), then TreeSHAP on the same
# samples. This is how Table S2 (12) ranks ASVs, here at genus level and on
# both platforms: a descriptive ranking of what the model relies on, with no
# split, so it is the same every time. Unmeasured cells stay blank. Numeric
# predictors are scaled, as in 09 and 12; this does not change the ranks.

import numpy as np
import pandas as pd
import shap
import xgboost as xgb
from scipy.stats import rankdata
from sklearn.preprocessing import StandardScaler

import figure_inputs as fi
from pipeline_settings import CAT_COLS, CORE_CHEM, OPS, TARGET, XGB_PARAMS, load_matrix, present

NAMES = {"PH": "pH", "TEMP": "Temperature", "HRT": "HRT",
         "Feed_Complexity": "Feed complexity", "Primary_Carbon_Signature": "Primary carbon signature",
         "g__Unclassified": "Unclassified (no genus call)"}


def display_name(feature):
    if feature in NAMES:
        return NAMES[feature]
    if feature.startswith("g__"):
        return feature[3:].replace("_", " ").strip()
    return feature


def fit():
    """Returns (X, shap_values, ranking). X holds the predictors as the model
    sees them before scaling (genera as relative abundance), shap_values the
    SHAP values in mM C caproate, and ranking every predictor's mean |SHAP| and
    rank (1 = largest; tied predictors share their average rank, as in 09)."""
    df, genus_cols = load_matrix(fi.matrix_path())
    features = sorted(set(present(OPS + CORE_CHEM + CAT_COLS, df) + genus_cols))
    X = df[features]
    Xs = X.copy()
    num = [c for c in features if c not in CAT_COLS]
    Xs[num] = StandardScaler().fit_transform(X[num])
    model = xgb.XGBRegressor(**XGB_PARAMS).fit(Xs, df[TARGET])
    shap_values = shap.TreeExplainer(model).shap_values(Xs)
    mean_abs = np.abs(shap_values).mean(axis=0)
    ranking = pd.DataFrame({"feature": features, "mean_abs_shap": mean_abs,
                            "rank": rankdata(-mean_abs, method="average")})
    ranking = ranking.sort_values(["rank", "feature"]).reset_index(drop=True)
    print(f"One XGBoost fit on all {len(df)} modelled samples ({df['Paper_ID'].nunique()} studies) with "
          f"{len(features)} predictors ({len(genus_cols)} genera); TreeSHAP on the same samples.")
    return X, shap_values, ranking
