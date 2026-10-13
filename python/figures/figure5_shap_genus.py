# figure5_shap_genus.py
# Figure 5: SHAP summary for the genera the model relies on most, from one
# XGBoost fit on every modelled sample of the corrected matrix with the
# headline's predictors and settings (all_samples_shap.py; the same fit as
# Figure S5). It shows the top 10 genera by mean |SHAP|, 10 as in 09's top-10
# criterion. g__Unclassified is left out, as in Figure 6: it is a bin of
# unassigned reads, not a taxon (its rank is printed). Run it with
# run_figures.py, which points it at the matrix.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata

import all_samples_shap
import figure_inputs as fi

plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = 9
plt.rcParams["axes.labelsize"] = 10
plt.rcParams["axes.titlesize"] = 10
plt.rcParams["xtick.labelsize"] = 8.5
plt.rcParams["ytick.labelsize"] = 8.5
plt.rcParams["legend.fontsize"] = 8.5
plt.rcParams["figure.titlesize"] = 11
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

N_TOP_GENERA = 10
NOT_A_TAXON = "g__Unclassified"


def generate_figure_5():
    """Figure 5 -- SHAP summary for the top genera, computed from a real
    XGBoost fit + shap.TreeExplainer on every modelled sample. No simulated
    SHAP or feature values anywhere in this function."""
    X, shap_values, ranking = all_samples_shap.fit()

    genus_ranked = ranking[ranking["feature"].str.startswith("g__")].reset_index(drop=True)
    genus_ranked["rank_among_genera"] = rankdata(-genus_ranked["mean_abs_shap"], method="average")
    genus_ranked = genus_ranked.rename(columns={"rank": "rank_among_all_predictors"})
    print("\n[Figure 5] Top 15 genera by mean |SHAP| (mM C):")
    print(genus_ranked.head(15).to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    top_genera = genus_ranked.loc[genus_ranked["feature"] != NOT_A_TAXON, "feature"].head(N_TOP_GENERA).tolist()
    print(f"[Figure 5] Shown: the top {N_TOP_GENERA} genera, {NOT_A_TAXON} left out: {', '.join(top_genera)}")

    feature_index = {feat: i for i, feat in enumerate(X.columns)}
    plot_rows = []
    for feat in top_genera:
        idx = feature_index[feat]
        feat_vals = X[feat].values.astype(float)
        sv = shap_values[:, idx]
        rng = feat_vals.max() - feat_vals.min()
        fv_norm = (feat_vals - feat_vals.min()) / rng if rng > 0 else np.zeros_like(feat_vals)
        clean_name = all_samples_shap.display_name(feat)
        for s, f in zip(sv, fv_norm):
            plot_rows.append({"Taxon": clean_name, "SHAP Value": s, "Feature Value": f})
    df_shap = pd.DataFrame(plot_rows)

    taxa_order = [all_samples_shap.display_name(t) for t in top_genera]
    fig, ax = plt.subplots(figsize=(7.2, 0.44 * len(taxa_order) + 1.8))
    y_positions = {t: i for i, t in enumerate(reversed(taxa_order))}

    rng = np.random.default_rng(42)  # jitter only, not data -- purely visual de-overlap
    sc = None
    for taxon in taxa_order:
        sub = df_shap[df_shap["Taxon"] == taxon]
        y_base = y_positions[taxon]
        y_jitter = y_base + rng.normal(0, 0.08, len(sub))
        sc = ax.scatter(sub["SHAP Value"], y_jitter, c=sub["Feature Value"], cmap="coolwarm",
                        s=16, alpha=0.8, edgecolors="none")

    ax.axvline(0, color="#666666", linestyle="--", linewidth=0.8)
    ax.set_yticks([y_positions[t] for t in taxa_order])
    ax.set_yticklabels(taxa_order, fontweight="bold", fontstyle="italic")
    ax.set_xlabel("SHAP Value (Impact on Predicted Caproate, mM C)", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    cbar = plt.colorbar(sc, ax=ax, aspect=20, pad=0.02)
    cbar.set_label("Relative Abundance (within-genus, 0-1)", fontweight="bold", rotation=270, labelpad=12)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(["Low", "High"])

    plt.title("Figure 5. SHAP Attribution of the Genera the Model Relies on Most", loc="left", fontweight="bold", pad=12)
    plt.tight_layout()
    fi.save(fig, "Figure_5_SHAP_Summary")
    plt.close()


if __name__ == "__main__":
    generate_figure_5()
