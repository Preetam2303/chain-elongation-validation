import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import xgboost as xgb
import shap
from sklearn.preprocessing import StandardScaler

MATRIX_PATH = "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"
OUTPUT_DIR = "../../figures_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def generate_figure_s5(max_display=15):
    """Figure S5 -- full SHAP ranking across chemistry, operational parameters,
    AND genera together (unlike main-text Figure 5, which is deliberately
    genus-only to match Section 3.6's taxonomic focus). This is the
    metabolite-inclusive view -- answers "how do chemistry and biology compare
    directly" using shap.summary_plot on a real global model fit."""
    df = pd.read_csv(MATRIX_PATH, low_memory=False)

    genus_cols = sorted([c for c in df.columns if c.startswith("g__")])
    df[genus_cols] = df[genus_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    row_sums = df[genus_cols].sum(axis=1).replace(0, 1)
    df[genus_cols] = df[genus_cols].div(row_sums, axis=0)

    TARGET = "Caproate"
    foundational_chems = ["Lactate", "Acetate", "Butyrate", "Ethanol", "Propionate"]
    op_params = ["PH", "TEMP", "HRT"]
    valid_features = [c for c in foundational_chems + op_params if c in df.columns] + genus_cols

    X = df[valid_features].apply(pd.to_numeric, errors="coerce").fillna(0)
    y = pd.to_numeric(df[TARGET], errors="coerce").fillna(0)

    def clean_label(col_name):
        return col_name.replace("g__", "").replace("_", " ") if col_name.startswith("g__") else col_name
    X.columns = [clean_label(c) for c in X.columns]

    scaler = StandardScaler()
    X_scaled = pd.DataFrame(scaler.fit_transform(X), columns=X.columns, index=X.index)

    xgb_params = {"n_estimators": 150, "learning_rate": 0.05, "max_depth": 3, "subsample": 0.8,
                  "colsample_bytree": 0.8, "random_state": 42, "n_jobs": -1, "tree_method": "hist"}
    model = xgb.XGBRegressor(**xgb_params).fit(X_scaled, y)
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_scaled)

    mean_shap = np.abs(shap_values).mean(axis=0)
    shap_df = pd.DataFrame({"Feature": X_scaled.columns, "SHAP": mean_shap}).sort_values("SHAP", ascending=False)
    n_chem_ops = len([c for c in shap_df.head(max_display)["Feature"] if c in foundational_chems + op_params])
    print(f"[Figure S5] Top {max_display} by real mean |SHAP| ({n_chem_ops} chemistry/ops, {max_display - n_chem_ops} genera):")
    print(shap_df.head(max_display).to_string(index=False))

    plt.figure(figsize=(8, 6.5))
    shap.summary_plot(shap_values, X_scaled, max_display=max_display, show=False)
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("SHAP value (impact on predicted caproate, standardized units)", fontsize=11, fontweight="bold")
    plt.title("Figure S5. Combined Chemistry, Operational, and Genus-Level SHAP Attribution",
              loc="left", fontweight="bold", fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_S5_Combined_SHAP.png"))
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_S5_Combined_SHAP.pdf"))
    plt.close()
    print(f"[OK] Figure S5 saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    generate_figure_s5()
