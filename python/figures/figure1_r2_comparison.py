# figure1_r2_comparison.py
# The master validation-collapse figure: every naive/ungrouped scheme and
# every properly grouped/zero-leakage scheme, one bar chart. Values are the
# same verified scalars already reported in Tables 3a and 3b -- a legitimate
# use of hardcoded data (these are confirmed summary statistics, not
# per-sample data that could be fabricated), the same pattern used for
# Figure 3 and Figure S1.
#
# Every value below is cross-checked against the manuscript tables:
#    0.985            Table 3a, single-study internal 80/20 (Duber et al. 2025)
#    0.915            Table 3a, single random 80/20 split (n=130, 4037 features)
#    0.835 +/- 0.080  Table 3a, shuffled 10-fold CV
#    0.602 +/- 0.273  Table 3a, random 5-fold CV (ASV-only, 4027 features)
#    0.046 +/- 0.158  Table 3b, GroupKFold by vessel, all reactors
#    0.036 +/- 0.197  Table 3b, GroupKFold by vessel, CSTR only
#   -0.095 +/- 0.387  Table 3b / Table 4 tier 4, LOSO full genus set
#   -0.170            Table 3b, platform transfer, relative abundance
#   -0.207            Table 3b, platform transfer, raw counts
#   -0.493 +/- 0.375  Table 3b, LOSO SHAP top-25 dynamic
#
# Error = 0.0 is used ONLY for schemes that are a single train/test split and
# therefore have no fold-to-fold SD to report.

import os
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import pandas as pd

# Relative to this script's own location, matching every other figure script.
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "..", "..", "figures_output")
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
    # Ordered by performance inside their groups, with errors mapped perfectly
    data = [
        {"Scheme": "Single-Study Internal 80/20 (Duber 2025)", "R2": 0.985, "Error": 0.0, "Type": "Naive / Ungrouped"},
        {"Scheme": "Single Random 80/20 Split", "R2": 0.915, "Error": 0.0, "Type": "Naive / Ungrouped"},
        {"Scheme": "Shuffled 10-Fold CV (Full Matrix)", "R2": 0.835, "Error": 0.080, "Type": "Naive / Ungrouped"},
        {"Scheme": "Random 5-Fold CV", "R2": 0.602, "Error": 0.273, "Type": "Naive / Ungrouped"},
        {"Scheme": "GroupKFold by Vessel (All Reactors)", "R2": 0.046, "Error": 0.158, "Type": "Grouped / Zero-Leakage"},
        {"Scheme": "GroupKFold by Vessel (CSTR Only)", "R2": 0.036, "Error": 0.197, "Type": "Grouped / Zero-Leakage"},
        {"Scheme": "LOSO (Full Genus Set, Tier 4)", "R2": -0.095, "Error": 0.387, "Type": "Grouped / Zero-Leakage"},
        {"Scheme": "Platform Transfer (Relative Abundance)", "R2": -0.170, "Error": 0.0, "Type": "Grouped / Zero-Leakage"},
        {"Scheme": "Platform Transfer (Raw Counts)", "R2": -0.207, "Error": 0.0, "Type": "Grouped / Zero-Leakage"},
        {"Scheme": "LOSO (SHAP Top-25 Dynamic)", "R2": -0.493, "Error": 0.375, "Type": "Grouped / Zero-Leakage"},
    ]

    # Reverse the dataframe so the highest R2 values plot at the top of the chart
    df = pd.DataFrame(data).iloc[::-1].reset_index(drop=True)

    # Widened figure slightly so the outside legend fits perfectly
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    colors = [CLR_NAIVE if t == "Naive / Ungrouped" else CLR_GROUPED for t in df["Type"]]

    bars = ax.barh(df["Scheme"], df["R2"], xerr=df["Error"], color=colors, height=0.68,
                   edgecolor="none", error_kw={'ecolor': '#333333', 'capsize': 3, 'elinewidth': 1})

    ax.axvline(0, color="#333333", linestyle="--", linewidth=0.9)
    ax.set_xlabel("Predictive Accuracy ($R^2$)", fontweight="bold")

    # WIDENED THE NEGATIVE LIMIT TO -1.45 TO PREVENT TEXT COLLISION ON THE Y-AXIS
    ax.set_xlim(-1.45, 1.25)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for i, bar in enumerate(bars):
        width = bar.get_width()
        err = df["Error"].iloc[i]

        # Shift text dynamically past the error bar
        x_pos = width + err + 0.03 if width >= 0 else width - err - 0.03

        ax.text(x_pos, bar.get_y() + bar.get_height() / 2, f"{width:.3f}",
                va="center", ha="left" if width >= 0 else "right",
                fontsize=7.5, fontweight="bold", color="#333333")

    # Added the Zero-Variance line back into the custom handles
    legend_elements = [
        Patch(facecolor=CLR_NAIVE, label="Naive / Ungrouped Schemes"),
        Patch(facecolor=CLR_GROUPED, label="Grouped / Zero-Leakage Schemes"),
        Line2D([0], [0], color="#333333", linestyle="--", linewidth=0.9, label="Zero-Variance ($R^2 = 0$)")
    ]

    # Pull the legend completely outside the right edge of the plot
    ax.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1.02, 0.5),
              frameon=True, edgecolor="#cccccc")

    plt.title("Figure 1. Predictive Performance Collapse Under Zero-Leakage Validation",
              loc="left", fontweight="bold", pad=12)

    # bbox_inches="tight" ensures the saved image expands to include the outside legend
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_1_Validation_Architectures.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_1_Validation_Architectures.pdf"), bbox_inches="tight")
    plt.close()
    print(f"[OK] Figure 1 saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    generate_figure_1()
