import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
from scipy.stats import spearmanr

MATRIX_PATH = "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"
OUTPUT_DIR = "../../figures_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = 9
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

CLR_NAIVE = "#E64B35"


def generate_figure_s4():
    """Figure S4 -- actual vs. predicted caproate for the Duber et al. (2025)
    B1 -> B2 transfer (Table 5; R2 = -5.869, Spearman rho = -0.569). Real
    train-only-scaled XGBoost fit and real predictions on the real held-out
    reactor -- the same leakage-controlled procedure used throughout this
    study (Section 2.6), and matching the exact feature-list and
    reactor-selection logic of the script that produced Table 5.

    Note on exact reproduction: XGBoost does not guarantee bit-identical
    output across library versions even with a fixed random_state -- internal
    tie-breaking and numerics can differ between releases. Run this in the
    same environment that produced Table 5 (i.e., on your machine, with your
    pinned package versions) to get -5.869 exactly; a different XGBoost
    version will give a close but not identical number, which is a known
    library characteristic, not a bug in this script."""
    df = pd.read_csv(MATRIX_PATH, low_memory=False)

    genus_cols = sorted([c for c in df.columns if c.startswith("g__")])
    df[genus_cols] = df[genus_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    row_sums = df[genus_cols].sum(axis=1).replace(0, 1)
    df[genus_cols] = df[genus_cols].div(row_sums, axis=0)

    TARGET = "Caproate"
    cat_cols = ["Feed_Complexity", "Primary_Carbon_Signature"]
    foundational_chems = ["Lactate", "Acetate", "Butyrate", "Ethanol", "Propionate"]
    op_params = ["PH", "TEMP", "HRT"]
    # sorted(set(...)) exactly, matching the script that produced the canonical
    # Table 5 numbers -- an unsorted list changes which columns XGBoost's
    # colsample_bytree samples at a given index even under a fixed seed
    valid_features = sorted(list(set(foundational_chems + op_params + cat_cols + genus_cols)))

    df_p = df[df["Paper_ID"] == "Duber_2024"].copy()  # Paper_ID label; reported in-text as Duber et al. 2025
    for c in ["Feed_Complexity", "Primary_Carbon_Signature"]:
        df_p[c] = df_p[c].astype("category")

    # exact reactor-assignment logic from the canonical script: whichever two
    # bioreactors appear first in row order become train/test, not a hardcoded
    # 'B1'/'B2' assumption. For Duber_2024 this resolves to B1 -> B2 either way
    # (confirmed), but matches the canonical source exactly for full fidelity.
    bioreactors = df_p["BIOREACTOR"].dropna().unique()
    br_train, br_test = bioreactors[0], bioreactors[1]
    train_df = df_p[df_p["BIOREACTOR"].astype(str).str.strip() == str(br_train).strip()].copy()
    test_df = df_p[df_p["BIOREACTOR"].astype(str).str.strip() == str(br_test).strip()].copy()

    X_train, X_test = train_df[valid_features].copy(), test_df[valid_features].copy()
    y_train = pd.to_numeric(train_df[TARGET], errors="coerce").fillna(0)
    y_test = pd.to_numeric(test_df[TARGET], errors="coerce").fillna(0)

    numeric_cols = [c for c in valid_features if c not in ["Feed_Complexity", "Primary_Carbon_Signature"]]
    for c in numeric_cols:
        X_train[c] = pd.to_numeric(X_train[c], errors="coerce").fillna(0)
        X_test[c] = pd.to_numeric(X_test[c], errors="coerce").fillna(0)

    scaler = StandardScaler()  # fit on train only -- no leakage
    X_train[numeric_cols] = scaler.fit_transform(X_train[numeric_cols])
    X_test[numeric_cols] = scaler.transform(X_test[numeric_cols])

    xgb_params = {"n_estimators": 150, "learning_rate": 0.05, "max_depth": 3, "subsample": 0.8,
                  "colsample_bytree": 0.8, "random_state": 42, "n_jobs": -1,
                  "enable_categorical": True, "tree_method": "hist"}
    model = xgb.XGBRegressor(**xgb_params).fit(X_train, y_train)
    y_pred = model.predict(X_test)

    r2 = r2_score(y_test, y_pred)
    rho, _ = spearmanr(y_test, y_pred)
    print(f"[Figure S4] Real recomputation: R2={r2:.3f}, Spearman rho={rho:.3f} (should match Table 5: -5.869, -0.569)")

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(y_test, y_pred, color=CLR_NAIVE, alpha=0.85, edgecolor="black", s=45, linewidth=0.6)

    lims = [min(y_test.min(), y_pred.min()) - 10, max(y_test.max(), y_pred.max()) + 10]
    ax.plot(lims, lims, "k--", linewidth=1.0, label="Perfect prediction ($y = x$)")

    ax.set_xlabel("Actual Caproate (mM C)", fontweight="bold")
    ax.set_ylabel("Predicted Caproate (mM C)", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper left", frameon=True, edgecolor="none", fontsize=8)
    ax.text(0.97, 0.03, f"$R^2 = {r2:.3f}$\nSpearman $\\rho = {rho:.3f}$", transform=ax.transAxes,
            fontsize=8.5, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#F8F9FA", edgecolor="#CCCCCC"))

    plt.title("Figure S4. Rank-Inverted Prediction (Duber et al. 2025, B1 \u2192 B2)",
              loc="left", fontweight="bold", pad=12, fontsize=9.5)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_S4_Duber_Scatter.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_S4_Duber_Scatter.pdf"), bbox_inches="tight")
    plt.close()
    print(f"[OK] Figure S4 saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    generate_figure_s4()
