# figure1_r2_comparison.py
# The master validation-collapse figure: every naive/ungrouped scheme and
# every properly grouped/zero-leakage scheme, one bar chart. Values are the
# same verified scalars already reported in Tables 3a and 3b -- a legitimate
# use of hardcoded data (these are confirmed summary statistics, not
# per-sample data that could be fabricated), the same pattern used for
# Figure 3 and Figure S1.

import os
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import pandas as pd

OUTPUT_DIR = "../../figures_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

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

CLR_NAIVE = "#E64B35"
CLR_GROUPED = "#3C5488"


def generate_figure_1():
    data = [
        {"Scheme": "Single-Study Internal 80/20 (Duber 2025)", "R2": 0.985, "Type": "Naive / Ungrouped"},
        {"Scheme": "Shuffled 10-Fold CV (Full Matrix)", "R2": 0.835, "Type": "Naive / Ungrouped"},
        {"Scheme": "Single Random 80/20 Split", "R2": 0.809, "Type": "Naive / Ungrouped"},
        {"Scheme": "Random 5-Fold CV", "R2": 0.668, "Type": "Naive / Ungrouped"},
        {"Scheme": "GroupKFold by Vessel (All Reactors)", "R2": 0.046, "Type": "Grouped / Zero-Leakage"},
        {"Scheme": "GroupKFold by Vessel (CSTR Only)", "R2": 0.036, "Type": "Grouped / Zero-Leakage"},
        {"Scheme": "LOSO (Full Genus Set, Tier 4)", "R2": -0.095, "Type": "Grouped / Zero-Leakage"},
        {"Scheme": "Platform Transfer (Relative Abundance)", "R2": -0.170, "Type": "Grouped / Zero-Leakage"},
        {"Scheme": "Platform Transfer (Raw Counts)", "R2": -0.207, "Type": "Grouped / Zero-Leakage"},
        {"Scheme": "LOSO (SHAP Top-25 Dynamic)", "R2": -0.493, "Type": "Grouped / Zero-Leakage"},
    ]

    df = pd.DataFrame(data)
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    colors = [CLR_NAIVE if t == "Naive / Ungrouped" else CLR_GROUPED for t in df["Type"]]
    bars = ax.barh(df["Scheme"], df["R2"], color=colors, height=0.68, edgecolor="none")

    ax.axvline(0, color="#333333", linestyle="--", linewidth=0.9, label="Zero-Variance ($R^2 = 0$)")
    ax.set_xlabel("Predictive Accuracy ($R^2$)", fontweight="bold")
    ax.set_xlim(-0.65, 1.1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for bar in bars:
        width = bar.get_width()
        x_pos = width + 0.02 if width >= 0 else width - 0.08
        ax.text(x_pos, bar.get_y() + bar.get_height() / 2, f"{width:.3f}",
                va="center", ha="left" if width >= 0 else "right",
                fontsize=7.5, fontweight="bold", color="#333333")

    legend_elements = [
        Patch(facecolor=CLR_NAIVE, label="Naive / Ungrouped Schemes"),
        Patch(facecolor=CLR_GROUPED, label="Grouped / Zero-Leakage Schemes"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", frameon=True, edgecolor="none")

    plt.title("Figure 1. Predictive Performance Collapse Under Zero-Leakage Validation",
              loc="left", fontweight="bold", pad=12)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_1_Validation_Architectures.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_1_Validation_Architectures.pdf"), bbox_inches="tight")
    plt.close()
    print(f"[OK] Figure 1 saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    generate_figure_1()
