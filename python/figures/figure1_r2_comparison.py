# figure1_r2_comparison.py
# The master validation-collapse figure: every naive/ungrouped scheme and
# every properly grouped/zero-leakage scheme, one bar chart. Every value is
# read from the corrected run's results tables (Tables 3a and 3b), so the
# figure cannot drift from them:
#   10_table3a_naive_splits.csv   Table 3a: the four naive schemes (seed 42)
#   04_groupkfold_by_vessel.csv   Table 3b: GroupKFold by vessel, all vessels and CSTRs only
#   02_loso_tables_3b_4.csv       Table 3b: LOSO with the core set (the headline) and
#                                 with the per-fold SHAP top-25 genera
#   06_platform_transfer.csv      Table 3b: platform transfer, relative abundance and raw counts
# Error bars are the SD across folds; a scheme that is a single train/test
# split has none. Run it with run_figures.py, which points it at those files.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import figure_inputs as fi

plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["mathtext.fontset"] = "custom"  # R2 and rho in the text font (Arial), not DejaVu Sans
plt.rcParams["mathtext.cal"] = "sans"  # unused; avoids a warning where no cursive font is installed
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
NAIVE, GROUPED = "Naive / Ungrouped", "Grouped / Zero-Leakage"


def scheme(label, kind, row, r2="r2_mean", sd="r2_sd"):
    err = row[sd] if sd in row.index else np.nan
    return {"Scheme": label, "R2": float(row[r2]), "Error": 0.0 if pd.isna(err) else float(err), "Type": kind}


def generate_figure_1():
    t3a = pd.read_csv(fi.results_file("10", "10_table3a_naive_splits.csv"))
    vessel = pd.read_csv(fi.results_file("04", "04_groupkfold_by_vessel.csv"))
    loso = pd.read_csv(fi.results_file("02", "02_loso_tables_3b_4.csv"))
    platform = pd.read_csv(fi.results_file("06", "06_platform_transfer.csv"))

    xgb_vessel = vessel[vessel["model"].str.startswith("XGBoost")]
    xgb_platform = platform[platform["model"].str.startswith("XGBoost")]
    data = [
        scheme("Single-Study Internal 80/20 (Duber 2025)", NAIVE,
               fi.one_row(t3a, scheme="Single-study internal 80/20 split (Duber 2024)")),
        scheme("Single Random 80/20 Split (ASVs)", NAIVE,
               fi.one_row(t3a, scheme="Single random 80/20 split (ASV level)")),
        scheme("Shuffled 10-Fold CV (Full Matrix)", NAIVE,
               fi.one_row(t3a, scheme="Shuffled 10-fold CV, per-fold SHAP top-25 genera")),
        scheme("Random 5-Fold CV (ASVs)", NAIVE,
               fi.one_row(t3a, scheme="Random 5-fold CV (ASV level)")),
        scheme("GroupKFold by Vessel (All Vessels)", GROUPED,
               fi.one_row(xgb_vessel, scheme="Mode 2: All Vessels (CSTR + Batch R1-R9)")),
        scheme("GroupKFold by Vessel (CSTRs Only)", GROUPED,
               fi.one_row(xgb_vessel, scheme="Mode 1: CSTR Only (B1, B2 Vessels)")),
        scheme("LOSO (Core Set; Headline)", GROUPED,
               fi.one_row(loso, column="full matrix",
                          model="Tier 3: + butyrate, propionate (core set; HEADLINE)")),
        scheme("LOSO (Per-Fold SHAP Top-25 Genera)", GROUPED,
               fi.one_row(loso, column="full matrix", model="LOSO, core set, per-fold SHAP top-25 genera")),
        scheme("Platform Transfer (Relative Abundance)", GROUPED,
               fi.one_row(xgb_platform, scheme="Platform transfer, relative abundance"), r2="r2"),
        scheme("Platform Transfer (Raw Counts)", GROUPED,
               fi.one_row(xgb_platform, scheme="Platform transfer, raw counts"), r2="r2"),
    ]
    # Ordered by R2 inside each group, naive schemes first
    df = pd.DataFrame(data)
    df = pd.concat([df[df["Type"] == NAIVE].sort_values("R2", ascending=False),
                    df[df["Type"] == GROUPED].sort_values("R2", ascending=False)])
    print("\n[Figure 1] Values plotted (R2 mean, SD across folds; 0 = single split):")
    print(df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))

    # Reverse the dataframe so the highest R2 values plot at the top of the chart
    df = df.iloc[::-1].reset_index(drop=True)

    # Widened figure slightly so the outside legend fits perfectly
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    colors = [CLR_NAIVE if t == NAIVE else CLR_GROUPED for t in df["Type"]]

    bars = ax.barh(df["Scheme"], df["R2"], xerr=df["Error"], color=colors, height=0.68,
                   edgecolor="none", error_kw={'ecolor': '#333333', 'capsize': 3, 'elinewidth': 1})

    ax.axvline(0, color="#333333", linestyle="--", linewidth=0.9)
    ax.set_xlabel("Predictive Accuracy ($R^2$)", fontweight="bold")

    low = (df["R2"] - df["Error"]).min()
    high = (df["R2"] + df["Error"]).max()
    ax.set_xlim(min(-1.0, low - 0.1), max(1.0, high + 0.1))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    labels = []
    for i, bar in enumerate(bars):
        width = bar.get_width()
        err = df["Error"].iloc[i]

        # Each value sits just past the end of its error bar (or bar), on the outer side
        end = width + err if width >= 0 else width - err
        labels.append(ax.annotate(f"{width:.3f}", (end, bar.get_y() + bar.get_height() / 2),
                                  xytext=(5 if width >= 0 else -5, 0), textcoords="offset points",
                                  va="center", ha="left" if width >= 0 else "right",
                                  fontsize=7.5, fontweight="bold", color="#333333"))

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

    # Widen the x-axis until every value label fits inside it. The labels have a
    # fixed size in points, so each widening makes them a little longer in R2
    # units; a few rounds settle it.
    for _ in range(4):
        fig.canvas.draw()
        boxes = [t.get_window_extent().transformed(ax.transData.inverted()) for t in labels]
        left, right = ax.get_xlim()
        ax.set_xlim(min(left, min(b.x0 for b in boxes) - 0.05), max(right, max(b.x1 for b in boxes) + 0.05))
    fi.save(fig, "Figure_1_Validation_Architectures")
    plt.close()


if __name__ == "__main__":
    generate_figure_1()
