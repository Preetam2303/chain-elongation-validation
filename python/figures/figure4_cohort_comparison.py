# figure4_cohort_comparison.py
# Figure 4: Illumina vs. Nanopore cohort comparison, built entirely from the
# corrected run's genus matrix (BIOTWIN_MATRIX, from 00b): the modelled samples
# only, so the Illumina cohort is the training cohort of the platform transfer
# (06). Unmeasured values are blank and left out, never counted as 0. Run it
# with run_figures.py, which points it at the matrix.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import figure_inputs as fi
from pipeline_settings import NANOPORE_ID, load_matrix

# ==============================================================================
# SAME JOURNAL STYLE AS FIGURES 1 & 3 -- kept identical for visual consistency
# ==============================================================================
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

CLR_ILLUMINA = "#3C5488"
CLR_NANOPORE = "#00A087"


def generate_figure_4():
    """Figure 4 -- Illumina vs. Nanopore cohort comparison, built entirely from
    the corrected matrix. No simulated values anywhere in this function."""
    df, genus_cols = load_matrix(fi.matrix_path())  # genera as relative abundance

    nanopore = df["Paper_ID"] == NANOPORE_ID
    ill_label, nano_label = f"Illumina\n(n={int((~nanopore).sum())})", f"Nanopore\n(n={int(nanopore.sum())})"
    cohort = np.where(nanopore, nano_label, ill_label)

    plot_df = pd.DataFrame({
        "Cohort": cohort,
        "pH": df["PH"],
        "Lactate (mM C)": df["Lactate"],
        "Acetate (mM C)": df["Acetate"],
        "Caproate (mM C)": df["Caproate"],
        "Caproiciproducens (%)": 100 * df["g__Caproiciproducens"],
        "Lacticaseibacillus (%)": 100 * df["g__Lacticaseibacillus"],
    })

    metrics = ["pH", "Lactate (mM C)", "Acetate (mM C)", "Caproate (mM C)",
               "Caproiciproducens (%)", "Lacticaseibacillus (%)"]

    fig, axes = plt.subplots(2, 3, figsize=(8.5, 5.5))
    axes = axes.flatten()
    palette = {ill_label: CLR_ILLUMINA, nano_label: CLR_NANOPORE}
    order = [ill_label, nano_label]

    np.random.seed(42)  # the strip plots' jitter only, fixed so every run draws the same figure
    for i, metric in enumerate(metrics):
        ax = axes[i]
        # hue=Cohort + legend=False avoids the palette-without-hue deprecation
        # warning that showed up in the original run's console output
        sns.boxplot(data=plot_df, x="Cohort", y=metric, hue="Cohort", order=order,
                    palette=palette, ax=ax, width=0.4, fliersize=0,
                    boxprops=dict(alpha=0.8), legend=False)
        sns.stripplot(data=plot_df, x="Cohort", y=metric, hue="Cohort", order=order,
                      palette=palette, ax=ax, size=3, jitter=0.2, alpha=0.5, legend=False)
        ax.set_xlabel("")
        ax.set_ylabel(metric, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.suptitle("Figure 4. Distributional & Compositional Shift Between Sequencing Platforms",
                 x=0.08, y=0.98, ha="left", fontweight="bold", fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fi.save(fig, "Figure_4_Cohort_Comparison")
    plt.close()

    # Printout for the text: per cohort, the samples with a value, mean and median
    summary = plot_df.groupby("Cohort")[metrics].agg(["count", "mean", "median"]).T.unstack()
    summary = summary.rename(columns=lambda c: c.replace("\n", " "), level=0)
    print("\n[Figure 4] Per cohort: samples with a value (count), mean and median:")
    print(summary.to_string(float_format=lambda x: f"{x:.2f}"))


if __name__ == "__main__":
    generate_figure_4()
