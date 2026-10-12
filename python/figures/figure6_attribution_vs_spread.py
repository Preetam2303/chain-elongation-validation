# figure6_attribution_vs_spread.py
# -----------------------------------------------------------------------------
# One point per genus. Horizontal: in how many of the study groups the genus
# is detected at all. Vertical: its median importance rank (mean |SHAP|)
# across 09's validation folds, rank 1 at the top.
#
# The message the figure carries:
#   - Top-right  = important AND found everywhere. A genuine universal driver
#                  would sit here.
#   - Top-left   = important but confined to one or two studies. These genera
#                  rank highly because they identify WHICH STUDY a sample came
#                  from, not because they drive caproate.
#   - Bottom     = rarely or never important.
# An asterisk marks a named genus that gets exactly zero attribution once the
# Nanopore study is held out (09's zero_when_no_nano): its importance comes
# from that study alone, wherever else it is detected.
#
# The flat row of points near rank ~130 is genera that received exactly zero
# attribution in most folds; tied genera share one averaged rank (script 09).
# g__Unclassified is excluded: it is a bin of unassigned reads, not a taxon.
# Named: every genus with median rank 10 or better, and Lacticaseibacillus,
# found in most samples but rarely important.
#
# Reads 09's guild_membership_shap_ranks.csv and its log (the fold count) from
# 09's run of record, and the corrected matrix (BIOTWIN_MATRIX). Run it with
# run_figures.py, which points it at those files.
# -----------------------------------------------------------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

import figure_inputs as fi
from pipeline_settings import load_matrix

plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = 9
plt.rcParams["axes.labelsize"] = 10
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

CLR_POS = "#00A087"    # abundance pushes predicted caproate up
CLR_NEG = "#E64B35"    # abundance pushes predicted caproate down
CLR_FLAT = "#B0B0B0"   # no consistent direction
DIR_THRESHOLD = 0.15   # |shap_dir| below this is treated as no consistent direction

LABEL_RANK = 10                            # name every genus with median rank this good or better
EXTRA_LABELS = ["g__Lacticaseibacillus"]
# (dx, dy) offset in points from the point, set for the run of record; others get (10, 0)
OFFSETS = {
    "g__Caproiciproducens":            (-16, -10),
    "g__Rummeliibacillus":             (10, 0),
    "g__Gudongella":                   (10, 2),
    "g__Microbacterium":               (10, 0),
    "g__Anaerobium":                   (10, 7),
    "g__Ochrobactrum":                 (10, -7),
    "g__Clostridium_sensu_stricto_12": (-14, 8),
    "g__Oscillibacter":                (-14, -9),
    "g__Sporolactobacillus":           (10, 7),
    "g__Monoglobus":                   (10, -8),
    "g__Lacticaseibacillus":           (10, 0),
}
SHORT_NAMES = {"g__Clostridium_sensu_stricto_12": "C. sensu stricto 12",
               "g__Clostridium_sensu_stricto_1": "C. sensu stricto 1"}


def main():
    ranks = pd.read_csv(fi.results_file("09", "guild_membership_shap_ranks.csv")).set_index("genus")
    n_folds, n_samples, n_genera = fi.log_number(
        "09", r"^GUILD MEMBERSHIP BY SHAP RANK\s+--\s+(\d+) folds on the (\d+)-sample, (\d+)-genus matrix$")
    df, genus_cols = load_matrix(fi.matrix_path())
    if (len(df), len(genus_cols)) != (n_samples, n_genera) or set(genus_cols) != set(ranks.index):
        raise SystemExit(f"09 ranked {n_genera:.0f} genera on {n_samples:.0f} samples; the matrix has "
                         f"{len(genus_cols)} genera and {len(df)} samples.")
    n_groups = df["Paper_ID"].nunique()

    rows = []
    for g in genus_cols:
        present = df[g] > 0
        if present.sum() == 0:
            continue
        if g == "g__Unclassified":     # catch-all bin of unassigned reads, not a taxon
            continue
        r = ranks.loc[g]
        rows.append({"genus": g,
                     "groups": df.loc[present, "Paper_ID"].nunique(),
                     "samples": int(present.sum()),
                     "median_rank": float(r["median_rank"]),
                     "shap_dir": float(r["shap_dir"]) if pd.notna(r["shap_dir"]) else 0.0,
                     "zero_when_no_nano": bool(r["zero_when_no_nano"])})
    d = pd.DataFrame(rows)
    d["named"] = (d["median_rank"] <= LABEL_RANK) | d["genus"].isin(EXTRA_LABELS)

    rng = np.random.default_rng(7)               # fixed jitter, reproducible
    d["x"] = d["groups"] + rng.uniform(-0.28, 0.28, len(d))
    # keep named points un-jittered so labels sit where the reader expects;
    # named genera at exactly the same place are spread a little sideways
    d.loc[d["named"], "x"] = d.loc[d["named"], "groups"].astype(float)
    for _, same in d[d["named"]].sort_values("genus").groupby(["groups", "median_rank"]):
        if len(same) > 1:
            d.loc[same.index, "x"] += np.linspace(-0.12, 0.12, len(same))

    def colour(v):
        if v > DIR_THRESHOLD:
            return CLR_POS
        if v < -DIR_THRESHOLD:
            return CLR_NEG
        return CLR_FLAT
    d["c"] = d["shap_dir"].apply(colour)
    d["s"] = 12 + 0.9 * d["samples"]

    fig, ax = plt.subplots(figsize=(8.8, 5.2))

    # interpretive zones: 'important' = median rank 5 or better; same threshold both zones
    ZONE = 5.5
    ax.add_patch(Rectangle((0.5, 0.8), 2.38, ZONE - 0.8, color="#F4E3DF", zorder=0, linewidth=0))
    ax.add_patch(Rectangle((n_groups - 0.62, 0.8), 1.12, ZONE - 0.8, color="#DDEFEA", zorder=0, linewidth=0))
    ax.text(0.62, 0.9, "Important, but found in only\none or two studies: identifies\nthe study, not the biology",
            fontsize=7.2, color="#8A3B2A", va="top", ha="left", style="italic")
    ax.text(n_groups + 0.42, 5.15, "Important AND\nfound everywhere",
            fontsize=7.2, color="#1E6B5A", va="bottom", ha="right", style="italic")

    unl = d[~d["named"]]
    ax.scatter(unl["x"], unl["median_rank"], s=unl["s"], c=unl["c"],
               alpha=0.45, edgecolors="none", zorder=2)
    lab = d[d["named"]]
    ax.scatter(lab["x"], lab["median_rank"], s=lab["s"], c=lab["c"],
               alpha=0.95, edgecolors="#222222", linewidths=0.7, zorder=3)

    for _, r in lab.iterrows():
        dx, dy = OFFSETS.get(r["genus"], (10, 0))
        text = SHORT_NAMES.get(r["genus"], r["genus"][3:].replace("_", " "))
        if r["zero_when_no_nano"]:
            text += "*"
        ha = "right" if dx < 0 else "left"
        ax.annotate(text, (r["x"], r["median_rank"]), xytext=(dx, dy),
                    textcoords="offset points", fontsize=8, style="italic",
                    ha=ha, va="center", zorder=4,
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="#555555",
                                    shrinkA=0, shrinkB=3))

    ax.set_yscale("log")
    ax.set_ylim(170, 0.8)                         # rank 1 at the top
    ax.set_yticks([1, 2, 5, 10, 25, 50, 100])
    ax.set_yticklabels(["1", "2", "5", "10", "25", "50", "100"])
    ax.set_xlim(0.5, n_groups + 0.5)
    ax.set_xticks(range(1, n_groups + 1))
    ax.set_xlabel(f"Number of study groups in which the genus is detected (of {n_groups})",
                  fontweight="bold")
    ax.set_ylabel(f"Median importance rank across {n_folds:.0f} validation folds\n(mean |SHAP|; 1 = most important)",
                  fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=CLR_POS, markersize=7,
               label="Raises predicted caproate"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=CLR_NEG, markersize=7,
               label="Lowers predicted caproate"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=CLR_FLAT, markersize=7,
               label="No consistent direction"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#999999", markersize=4,
               label="Size = samples detected"),
        Line2D([0], [0], linestyle="none", label="* No attribution once the\n   Nanopore study is held out"),
    ]
    ax.legend(handles=legend, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              fontsize=7.2, frameon=True, edgecolor="#cccccc", framealpha=0.95)

    plt.title("Figure 6. Model importance tracks study identity for every genus\n"
              "except the known chain elongator",
              loc="left", fontweight="bold", fontsize=10, pad=10)
    plt.tight_layout()
    fi.save(fig, "Figure_6_Attribution_vs_Spread")
    plt.close()

    print(f"\n[Figure 6] {len(d)} genera plotted ({n_folds:.0f} folds, {n_groups} study groups). Named genera:")
    print(lab[["genus", "groups", "samples", "median_rank", "shap_dir", "zero_when_no_nano"]]
          .sort_values("median_rank").to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    hi = d[d["median_rank"] <= 10]
    lo = d[d["median_rank"] > 25]
    print(f"     genera with median rank <=10: mean groups present = {hi['groups'].mean():.2f} (n={len(hi)})")
    print(f"     genera with median rank  >25: mean groups present = {lo['groups'].mean():.2f} (n={len(lo)})")


if __name__ == "__main__":
    main()
