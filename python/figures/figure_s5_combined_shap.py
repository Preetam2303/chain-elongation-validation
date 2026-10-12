# figure_s5_combined_shap.py
# Figure S5: full SHAP ranking across chemistry, operational parameters,
# the categorical descriptors AND genera together (unlike main-text Figure 5,
# which is genus-only). It is the same fit as Figure 5 (all_samples_shap.py):
# one XGBoost fit on every modelled sample with the headline's predictors and
# settings. Colours show each predictor's measured value (shap's own
# beeswarm); the categorical descriptors and unmeasured cells are grey. Run it
# with run_figures.py, which points it at the matrix.

import matplotlib.pyplot as plt
import numpy as np
import shap

import all_samples_shap
import figure_inputs as fi

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
    """Figure S5 -- the metabolite-inclusive view: how chemistry and biology
    compare directly in the same real global model fit as Figure 5."""
    X, shap_values, ranking = all_samples_shap.fit()
    top = ranking.head(max_display)
    n_genera = int(top["feature"].str.startswith("g__").sum())
    print(f"\n[Figure S5] Top {max_display} predictors by mean |SHAP| (mM C) "
          f"({max_display - n_genera} chemistry, operational or categorical; {n_genera} genera):")
    print(top.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    X_shown = X.rename(columns=all_samples_shap.display_name)
    plt.figure(figsize=(8, 6.5))
    shap.summary_plot(shap_values, X_shown, max_display=max_display, show=False,
                      rng=np.random.default_rng(42))  # drawing order of the dots only
    fig, ax = plt.gcf(), plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel("SHAP value (impact on predicted caproate, mM C)", fontsize=11, fontweight="bold")
    plt.title("Figure S5. Combined Chemistry, Operational, and\nGenus-Level SHAP Attribution",
              loc="left", fontweight="bold", fontsize=11)
    plt.tight_layout()
    fi.save(fig, "Figure_S5_Combined_SHAP")
    plt.close()


if __name__ == "__main__":
    generate_figure_s5()
