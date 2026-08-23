import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ==============================================================================
# CONFIG -- adjust MATRIX_PATH / OUTPUT_DIR to your local paths
# ==============================================================================
MATRIX_PATH = "../../data/BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv"
OUTPUT_DIR = "../../figures_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

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
    the real 184-sample matrix. No simulated values anywhere in this function."""
    df = pd.read_csv(MATRIX_PATH, low_memory=False)

    # relative abundance, same TSS procedure used throughout the pipeline
    genus_cols = sorted([c for c in df.columns if c.startswith("g__")])
    df[genus_cols] = df[genus_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    row_sums = df[genus_cols].sum(axis=1).replace(0, 1)
    rel = df[genus_cols].div(row_sums, axis=0) * 100  # as percent

    cohort = np.where(df["Paper_ID"] == "Hanna_2025", "Nanopore (n=54)", "Illumina (n=130)")

    plot_df = pd.DataFrame({
        "Cohort": cohort,
        "pH": pd.to_numeric(df["PH"], errors="coerce"),
        "Lactate (mM C)": pd.to_numeric(df["Lactate"], errors="coerce"),
        "Acetate (mM C)": pd.to_numeric(df["Acetate"], errors="coerce"),
        "Caproate (mM C)": pd.to_numeric(df["Caproate"], errors="coerce"),
        "Caproiciproducens (%)": rel["g__Caproiciproducens"],
        "Lactobacillus (%)": rel["g__Lactobacillus"],
        "Lacticaseibacillus (%)": rel["g__Lacticaseibacillus"],
    })

    metrics = ["pH", "Lactate (mM C)", "Acetate (mM C)", "Caproate (mM C)",
               "Caproiciproducens (%)", "Lactobacillus (%)", "Lacticaseibacillus (%)"]

    fig, axes = plt.subplots(2, 4, figsize=(11, 5.5))
    axes = axes.flatten()
    palette = {"Illumina (n=130)": CLR_ILLUMINA, "Nanopore (n=54)": CLR_NANOPORE}
    order = ["Illumina (n=130)", "Nanopore (n=54)"]

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

    fig.delaxes(axes[7])
    plt.suptitle("Figure 4. Distributional & Compositional Shift Between Sequencing Platforms",
                 x=0.08, y=0.98, ha="left", fontweight="bold", fontsize=11)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_4_Cohort_Comparison.png"), dpi=300)
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_4_Cohort_Comparison.pdf"))
    plt.close()

    # sanity-check printout so the console confirms these are the real, reported numbers
    print("[Figure 4] Real per-cohort means (should match manuscript Section 3.5):")
    print(plot_df.groupby("Cohort")[metrics].mean().round(2).to_string())
    print(f"[OK] Figure 4 saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    generate_figure_4()
