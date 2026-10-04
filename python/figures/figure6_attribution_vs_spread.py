# figure6_attribution_vs_spread.py
# -----------------------------------------------------------------------------
# One point per genus. Horizontal: in how many of the seven study groups the
# genus is detected at all. Vertical: its median importance rank (mean |SHAP|)
# across 22 validation folds, rank 1 at the top.
#
# The message the figure carries:
#   - Top-right  = important AND found everywhere. A genuine universal driver
#                  would sit here. Only Caproiciproducens does.
#   - Top-left   = important but confined to one or two studies. These genera
#                  rank highly because they identify WHICH STUDY a sample came
#                  from, not because they drive caproate. Several receive
#                  exactly zero attribution once their own cohort is removed.
#   - Bottom     = rarely or never important.
#
# The flat row of points near rank ~105 is genera that received exactly zero
# attribution in most folds; tied genera share one averaged rank (script 09).
# g__Unclassified is excluded: it is a bin of unassigned reads, not a taxon.
#
# Reads results/guild_membership_shap_ranks.csv (from script 09) and the final
# data matrix. Run from inside python/figures/.
# -----------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "data", "BIOTWIN_FINAL_GRAND_MERGE_Substrates.csv")
RANKS = os.path.join(HERE, "..", "..", "results", "guild_membership_shap_ranks.csv")
OUTPUT_DIR = os.path.join(HERE, "..", "..", "figures_output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

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

# label text, and (dx, dy) offset in points from the point
LABELS = {
    "g__Caproiciproducens":            ("Caproiciproducens",        (-16, -10)),
    "g__Oscillibacter":                ("Oscillibacter",            (-12, -14)),
    "g__Clostridium_sensu_stricto_12": ("C. sensu stricto 12",      (-14, 10)),
    "g__Clostridium_sensu_stricto_1":  ("C. sensu stricto 1",       (10, -2)),
    "g__Gudongella":                   ("Gudongella",               (10, 4)),
    "g__Acinetobacter":                ("Acinetobacter",            (10, 7)),
    "g__Anaerobium":                   ("Anaerobium",               (10, -2)),
    "g__Lacticaseibacillus":           ("Lacticaseibacillus",       (10, -3)),
}


def main():
    df = pd.read_csv(DATA)
    df["Paper_ID"] = df["Paper_ID"].astype(str).str.strip()
    genus_cols = [c for c in df.columns if str(c).startswith("g__")]
    df[genus_cols] = df[genus_cols].fillna(0)

    ranks = pd.read_csv(RANKS).set_index("genus")
    rows = []
    for g in genus_cols:
        present = df[g] > 0
        if present.sum() == 0 or g not in ranks.index:
            continue
        if g == "g__Unclassified":     # catch-all bin of unassigned reads, not a taxon
            continue
        r = ranks.loc[g]
        rows.append({"genus": g,
                     "groups": df.loc[present, "Paper_ID"].nunique(),
                     "samples": int(present.sum()),
                     "median_rank": float(r["median_rank"]),
                     "shap_dir": float(r["shap_dir"]) if pd.notna(r["shap_dir"]) else 0.0})
    d = pd.DataFrame(rows)

    rng = np.random.default_rng(7)               # fixed jitter, reproducible
    d["x"] = d["groups"] + rng.uniform(-0.28, 0.28, len(d))
    # keep labelled points un-jittered so labels sit where the reader expects
    d.loc[d["genus"].isin(LABELS), "x"] = d.loc[d["genus"].isin(LABELS), "groups"]

    def colour(v):
        if v > DIR_THRESHOLD:
            return CLR_POS
        if v < -DIR_THRESHOLD:
            return CLR_NEG
        return CLR_FLAT
    d["c"] = d["shap_dir"].apply(colour)
    d["s"] = 12 + 0.9 * d["samples"]

    fig, ax = plt.subplots(figsize=(7.2, 5.2))

    # interpretive zones
    ZONE = 5.5   # 'important' = median rank 5 or better; same threshold both zones
    ax.axhspan(0.8, ZONE, xmin=0.0, xmax=0.34, color="#F4E3DF", zorder=0)
    ax.axhspan(0.8, ZONE, xmin=0.84, xmax=1.0, color="#DDEFEA", zorder=0)
    ax.text(0.62, 5.15, "Important, but found in only\none or two studies: identifies\nthe study, not the biology",
            fontsize=7.2, color="#8A3B2A", va="bottom", ha="left", style="italic")
    ax.text(7.42, 5.15, "Important AND\nfound everywhere",
            fontsize=7.2, color="#1E6B5A", va="bottom", ha="right", style="italic")

    unl = d[~d["genus"].isin(LABELS)]
    ax.scatter(unl["x"], unl["median_rank"], s=unl["s"], c=unl["c"],
               alpha=0.45, edgecolors="none", zorder=2)
    lab = d[d["genus"].isin(LABELS)]
    ax.scatter(lab["x"], lab["median_rank"], s=lab["s"], c=lab["c"],
               alpha=0.95, edgecolors="#222222", linewidths=0.7, zorder=3)

    for _, r in lab.iterrows():
        text, (dx, dy) = LABELS[r["genus"]]
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
    ax.set_xlim(0.5, 7.5)
    ax.set_xticks(range(1, 8))
    ax.set_xlabel("Number of study groups in which the genus is detected (of 7)",
                  fontweight="bold")
    ax.set_ylabel("Median importance rank across 22 validation folds\n(mean |SHAP|; 1 = most important)",
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
    ]
    ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.555, 1.0),
              fontsize=7.2, frameon=True, edgecolor="#cccccc", framealpha=0.95)

    plt.title("Figure 6. Model importance tracks study identity for every genus\n"
              "except the known chain elongator",
              loc="left", fontweight="bold", fontsize=10, pad=10)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(os.path.join(OUTPUT_DIR, f"Figure_6_Attribution_vs_Spread.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close()

    hi = d[d["median_rank"] <= 10]
    lo = d[d["median_rank"] > 25]
    print(f"[OK] Figure 6 saved to {OUTPUT_DIR}")
    print(f"     genera with median rank <=10: mean groups present = {hi['groups'].mean():.2f} (n={len(hi)})")
    print(f"     genera with median rank  >25: mean groups present = {lo['groups'].mean():.2f} (n={len(lo)})")


if __name__ == "__main__":
    main()
