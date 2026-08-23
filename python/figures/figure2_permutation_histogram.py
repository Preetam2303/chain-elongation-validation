import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

OUTPUT_DIR = "../../figures_output"
RESULTS_DIR = "../../results"
PERMUTATION_RESULTS_CSV = os.path.join(RESULTS_DIR, "permutation_null_distribution_999.csv")
TRUE_OBSERVED_R2 = -0.095

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


def generate_figure_2(clip_at=-2.2):
    """Figure 2 -- label-permutation null distribution, loaded directly from the
    CSV saved by permutation_corrected.py. Every value plotted is a real,
    independently-computed LOSO result from a shuffled target; nothing here is
    simulated. The p-value is recomputed from the loaded data rather than
    hardcoded, so the figure can never drift out of sync with the numbers.

    Your real 999-run distribution has one genuine catastrophic outlier
    (R2 = -50.9) plus a handful of others below -2 -- the same kind of
    single-fold collapse seen elsewhere in this project (e.g. Duber 2025's
    -5.869). Left unclipped, matplotlib sizes the x-axis to fit that one point,
    which crushes the other ~990 values into a sliver -- that's the "one bar"
    you saw. The fix is to clip the DISPLAY range and say so on the plot, not
    to drop the point: the p-value below is still computed on the full,
    unclipped data.
    """
    perm_df = pd.read_csv(PERMUTATION_RESULTS_CSV)
    perm_r2_full = perm_df["permuted_r2"].values
    n_total = len(perm_r2_full)

    # statistics computed on the REAL, unclipped data
    k_beats_or_ties = int((perm_r2_full >= TRUE_OBSERVED_R2).sum())
    p_value = (k_beats_or_ties + 1) / (n_total + 1)
    n_clipped = int((perm_r2_full < clip_at).sum())
    min_val = perm_r2_full.min()

    # clip ONLY for display
    perm_r2_display = np.clip(perm_r2_full, a_min=clip_at, a_max=None)

    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    sns.histplot(perm_r2_display, bins=30, kde=False, color="#8491B4", edgecolor="white",
                 ax=ax, alpha=0.8)

    ax.axvline(TRUE_OBSERVED_R2, color=CLR_NAIVE, linestyle="--", linewidth=2.0,
               label=f"Observed LOSO $R^2$ ({TRUE_OBSERVED_R2:.3f})")
    ax.text(
        0.55, 0.75,
        f"Observed $R^2$: {TRUE_OBSERVED_R2:.3f}\nPermutations: {n_total}\n"
        f"Empirical $p = {p_value:.3f}$\n({n_total - k_beats_or_ties}/{n_total} permutations < true)",
        transform=ax.transAxes, fontsize=8, verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#F8F9FA", edgecolor="#CCCCCC", alpha=0.9),
    )

    # x-axis tick labels: show "<= clip_at" for the leftmost bin instead of a
    # misleading exact number, and note how many real points sit past it
    xticks = ax.get_xticks()
    xticklabels = [f"\u2264{clip_at:g}" if t <= clip_at else f"{t:g}" for t in xticks]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.set_xlim(clip_at - 0.15, max(0.3, perm_r2_full.max() + 0.1))

    ax.annotate(
        f"{n_clipped} permutations < {clip_at:g}\n(min = {min_val:.1f})",
        xy=(clip_at, ax.get_ylim()[1] * 0.9), xytext=(clip_at + 0.15, ax.get_ylim()[1] * 0.9),
        fontsize=7, color="#555555", ha="left", va="top",
        arrowprops=dict(arrowstyle="->", color="#555555", lw=0.8),
    )

    ax.set_xlabel("Permuted $R^2$ Score (display clipped at " + f"{clip_at:g}" + "; statistics use full data)", fontweight="bold", fontsize=8.5)
    ax.set_ylabel("Frequency (Number of Shuffles)", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper left", frameon=True, edgecolor="none")

    plt.title(f"Figure 2. Label-Permutation Null Distribution ({n_total} Iterations)",
              loc="left", fontweight="bold", pad=12)
    plt.tight_layout()

    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_2_Permutation_Control.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_2_Permutation_Control.pdf"), bbox_inches="tight")
    plt.close()

    print(f"[Figure 2] Loaded {n_total} real permuted R2 values from {PERMUTATION_RESULTS_CSV}")
    print(f"[Figure 2] Full-data mean={perm_r2_full.mean():.3f}, SD={perm_r2_full.std():.3f}, p={p_value:.4f}")
    print(f"[Figure 2] {n_clipped} points below display floor of {clip_at:g} (min={min_val:.2f}); shown via annotation, not hidden")
    print(f"[OK] Figure 2 saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    generate_figure_2()
