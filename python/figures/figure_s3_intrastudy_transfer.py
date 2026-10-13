# figure_s3_intrastudy_transfer.py
# Figure S3: Table 5 as a bar chart ordered by R2, the reactor-pair transfers
# within each study, R2 next to the Spearman rho between measured and predicted
# caproate (on sampling-day means, so same-day replicates count once). Every
# value is read from 05's run of record (05_table5_reactor_pairs.csv). Run it
# with run_figures.py, which points it at that file.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figure_inputs as fi

plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["mathtext.fontset"] = "custom"  # R2 and rho in the text font (Arial), not DejaVu Sans
plt.rcParams["mathtext.cal"] = "sans"  # unused; avoids a warning where no cursive font is installed
plt.rcParams["font.size"] = 9
plt.rcParams["axes.labelsize"] = 10
plt.rcParams["axes.titlesize"] = 10
plt.rcParams["xtick.labelsize"] = 8
plt.rcParams["ytick.labelsize"] = 8.5
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

CLR_NAIVE = "#E64B35"
CLR_GROUPED = "#3C5488"


def generate_figure_s3():
    t5 = pd.read_csv(fi.results_file("05", "05_table5_reactor_pairs.csv"))
    t5 = t5.sort_values("R2", ascending=False).reset_index(drop=True)
    print("\n[Figure S3] Values plotted (Table 5, ordered by R2):")
    print(t5.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    studies = [f"{fi.STUDY_LABELS[p]}\n({a} \u2192 {b})"
               for p, a, b in zip(t5["Paper"], t5["Train_Reactor"], t5["Test_Reactor"])]
    r2_vals = t5["R2"].tolist()
    rho_vals = t5["Spearman_rho_unit_means"].tolist()

    x = np.arange(len(studies))
    width = 0.35

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    bars_r2 = ax.bar(x - width / 2, r2_vals, width, label="$R^2$",
                     color=[CLR_GROUPED if v >= 0 else CLR_NAIVE for v in r2_vals],
                     edgecolor="black", linewidth=0.6)
    bars_rho = ax.bar(x + width / 2, rho_vals, width, label="Spearman $\\rho$ (sampling-day means)",
                      color="#8491B4", edgecolor="black", linewidth=0.6, hatch="//")

    ax.axhline(0, color="#333333", linewidth=1.0)
    ax.set_ylabel("Value (transfer $R^2$ or rank correlation)", fontweight="bold", labelpad=10)
    ax.set_ylim(min(min(r2_vals), min(rho_vals)) - 1.1, 1.5)
    ax.set_xticks(x)
    ax.set_xticklabels(studies)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower left", frameon=True, edgecolor="none")

    for bar in list(bars_r2) + list(bars_rho):
        yval = bar.get_height()
        offset = 0.2 if yval >= 0 else -0.35
        ax.text(bar.get_x() + bar.get_width() / 2, yval + offset, f"{yval:.3f}",
                ha="center", va="bottom" if yval >= 0 else "top", fontweight="bold", fontsize=7.5)

    plt.title("Figure S3 (optional). Intra-Study Reactor-Pair Transferability: $R^2$ vs. Rank Correlation (Table 5)",
              loc="left", fontweight="bold", pad=12, fontsize=9)
    plt.tight_layout()
    fi.save(fig, "Figure_S3_Intrastudy_Transfer")
    plt.close()


if __name__ == "__main__":
    generate_figure_s3()
