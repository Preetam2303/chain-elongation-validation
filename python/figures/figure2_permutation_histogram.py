# figure2_permutation_histogram.py
# Figure 2: the label-permutation controls of the headline LOSO (03). Two
# shuffles of caproate by sampling day, 999 times each: across all studies,
# and within each study. Each is shown by mean R2 (higher is better) and by
# mean RMSE (lower is better), in four panels. Every value comes from 03's run
# of record:
#   permutation_null_distribution_999.csv   every shuffle's result, across studies
#   permutation_null_within_study_999.csv   every shuffle's result, within each study
#   03_permutation_tests.csv                observed values, k and p
# k and p are recomputed here from the shuffles with 03's rule (a shuffle
# counts when its R2 reaches the observed R2, rounded to 3 decimals as in 03,
# or when its RMSE is as low as the observed one or lower; p = (k + 1) / 1000),
# and the script stops if they differ from 03's table.
#
# A few shuffles give extreme R2 values (single-fold collapses). Left
# unclipped, the axis would be sized to fit them and crush the other shuffles
# into a sliver, so the display is clipped at the 1st percentile (rounded down
# to 0.5) and the panel says how many shuffles lie below it. The statistics use
# every shuffle, unclipped.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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
CLR_NULL = "#8491B4"
CLR_BENCHMARK = "#555555"

SHUFFLES = [
    ("across studies", "permutation_null_distribution_999.csv", "Caproate shuffled by sampling day across studies"),
    ("within each study", "permutation_null_within_study_999.csv", "Caproate shuffled by sampling day within each study"),
]


def check(name, ours, theirs):
    if not np.isclose(ours, theirs, rtol=0, atol=1e-9):
        raise SystemExit(f"{name}: {ours} recomputed here, {theirs} in 03_permutation_tests.csv.")


def panel(ax, null, observed, k, p, metric, benchmark=None):
    n = len(null)
    if metric == "r2":
        clip_at = np.floor(np.percentile(null, 1) * 2) / 2
        shown = np.clip(null, clip_at, None)
        n_clipped = int((null < clip_at).sum())
    else:
        clip_at, shown, n_clipped = None, null, 0
    sns.histplot(shown, bins=30, kde=False, color=CLR_NULL, edgecolor="white", ax=ax, alpha=0.8)
    ax.set_ylim(0, 1.5 * max(bar.get_height() for bar in ax.patches))  # room for the box above the bars
    ax.axvline(observed, color=CLR_NAIVE, linestyle="--", linewidth=2.0)
    if benchmark is not None:
        ax.axvline(benchmark, color=CLR_BENCHMARK, linestyle=":", linewidth=1.5)

    if metric == "r2":
        obs_text, better = f"Observed $R^2$: {observed:.3f}", "as high or higher"
    else:
        obs_text, better = f"Observed RMSE: {observed:.1f} mM C", "as low or lower"
    lines = [obs_text, f"{k} of {n} shuffles {better}", f"Empirical $p = {p:.3f}$"]
    if n_clipped:
        lines += [f"{n_clipped} shuffles below {clip_at:g} drawn at {clip_at:g}", f"(lowest {null.min():.1f})"]

    lo = min(shown.min(), observed) if benchmark is None else min(shown.min(), observed, benchmark)
    hi = max(shown.max(), observed) if benchmark is None else max(shown.max(), observed, benchmark)
    pad = 0.05 * (hi - lo)
    ax.set_xlim(lo - pad, hi + pad)
    # the box goes on the side away from the observed line
    left = observed > (lo + hi) / 2
    ax.text(0.03 if left else 0.97, 0.95, "\n".join(lines),
            transform=ax.transAxes, fontsize=7.5, va="top", ha="left" if left else "right",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#F8F9FA", edgecolor="#CCCCCC", alpha=0.9))
    if clip_at is not None and n_clipped:
        ticks = [t for t in ax.get_xticks() if lo - pad <= t <= hi + pad]
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"\u2264{clip_at:g}" if t <= clip_at else f"{t:g}" for t in ticks])
    ax.set_ylabel("Number of shuffles", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return n_clipped, clip_at


def generate_figure_2():
    tests = pd.read_csv(fi.results_file("03", "03_permutation_tests.csv"))
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.4))
    print("\n[Figure 2] Values plotted (recomputed from every shuffle; checked against 03_permutation_tests.csv):")
    for row_axes, (shuffle, null_file, title) in zip(axes, SHUFFLES):
        null = pd.read_csv(fi.results_file("03", null_file))
        t = fi.one_row(tests, shuffle=shuffle)
        n = len(null)
        check(f"{shuffle}: number of shuffles", n, t["n_shuffles"])
        r2, rmse = null["permuted_r2"].to_numpy(), null["permuted_rmse"].to_numpy()
        k_r2 = int((r2 >= t["observed_r2"]).sum())
        k_rmse = int((rmse <= t["observed_rmse"]).sum())
        p_r2, p_rmse = (k_r2 + 1) / (n + 1), (k_rmse + 1) / (n + 1)
        for name, ours, theirs in [("k (R2)", k_r2, t["k_r2"]), ("p (R2)", p_r2, t["p_r2"]),
                                   ("k (RMSE)", k_rmse, t["k_rmse"]), ("p (RMSE)", p_rmse, t["p_rmse"]),
                                   ("null mean R2", r2.mean(), t["null_r2_mean"]),
                                   ("null SD R2", r2.std(), t["null_r2_sd"]),
                                   ("null mean RMSE", rmse.mean(), t["null_rmse_mean"]),
                                   ("null SD RMSE", rmse.std(), t["null_rmse_sd"])]:
            check(f"{shuffle}: {name}", ours, theirs)

        n_clipped, clip_at = panel(row_axes[0], r2, t["observed_r2"], k_r2, p_r2, "r2")
        panel(row_axes[1], rmse, t["observed_rmse"], k_rmse, p_rmse, "rmse",
              benchmark=t["training_mean_benchmark_rmse"])
        row_axes[0].set_title(title, loc="left", fontsize=9, fontweight="bold")
        print(f"  {shuffle}: {n} shuffles")
        print(f"    R2:   observed {t['observed_r2']:.3f} | null mean {r2.mean():.3f} (SD {r2.std():.3f}), "
              f"median {np.median(r2):.3f} | k = {k_r2}, p = {p_r2:.3f} | "
              f"display clipped at {clip_at:g} ({n_clipped} below, min {r2.min():.2f})")
        print(f"    RMSE: observed {t['observed_rmse']:.1f} mM C | null mean {rmse.mean():.1f} (SD {rmse.std():.1f}) "
              f"| k = {k_rmse}, p = {p_rmse:.3f} | training-mean benchmark {t['training_mean_benchmark_rmse']:.1f}")

    axes[1][0].set_xlabel("Mean LOSO $R^2$ of the shuffle (display clipped; statistics use every shuffle)",
                          fontweight="bold", fontsize=8.5)
    axes[1][1].set_xlabel("Mean LOSO RMSE of the shuffle (mM C)", fontweight="bold", fontsize=8.5)
    handles = [plt.Line2D([0], [0], color=CLR_NAIVE, linestyle="--", linewidth=2.0, label="Observed (true labels)"),
               plt.Line2D([0], [0], color=CLR_BENCHMARK, linestyle=":", linewidth=1.5,
                          label="Training-mean benchmark (RMSE)")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(f"Figure 2. Label-Permutation Null Distributions ({n} Shuffles Each)",
                 x=0.02, ha="left", fontweight="bold", fontsize=11)
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    fi.save(fig, "Figure_2_Permutation_Control")
    plt.close()


if __name__ == "__main__":
    generate_figure_2()
