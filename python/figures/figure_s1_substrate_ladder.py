# figure_s1_substrate_ladder.py
# Figure S1: the substrate ladder of Table 4, LOSO R2 per tier for the
# Illumina-only studies and for the full matrix (both platforms), with the SD
# across held-out studies. Every value is read from 02's run of record
# (02_loso_tables_3b_4.csv). The full matrix stops at tier 3, the headline core
# set: the tier-4 compounds were measured only in the Illumina studies. Run it
# with run_figures.py, which points it at that file.

import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import figure_inputs as fi

plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = 9
plt.rcParams["axes.labelsize"] = 10
plt.rcParams["axes.titlesize"] = 10
plt.rcParams["xtick.labelsize"] = 8.5
plt.rcParams["ytick.labelsize"] = 8.5
plt.rcParams["legend.fontsize"] = 8.5
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

CLR_GROUPED = "#3C5488"
CLR_ACCENT = "#00A087"
LABELS = {
    1: "Tier 1:\nOps + Genera",
    2: "Tier 2:\n+ Foundational Feeds",
    3: "Tier 3:\n+ Intermediates",
    4: "Tier 4:\n+ Downstream Proxies",
}


def tiers(table, column):
    rows = table[(table["column"] == column) & table["model"].str.match(r"Tier \d:")].copy()
    rows["tier"] = rows["model"].map(lambda m: int(re.match(r"Tier (\d):", m).group(1)))
    if rows["tier"].duplicated().any():
        raise SystemExit(f"More than one row per tier for {column} in 02_loso_tables_3b_4.csv.")
    return rows.set_index("tier").sort_index()


def generate_figure_s1():
    table = pd.read_csv(fi.results_file("02", "02_loso_tables_3b_4.csv"))
    illumina, full = tiers(table, "Illumina-only"), tiers(table, "full matrix")
    print("\n[Figure S1] Values plotted (LOSO R2 mean, SD across held-out studies):")
    for name, t in [("Illumina-only", illumina), ("Full matrix", full)]:
        for tier, r in t.iterrows():
            print(f"  {name:13s} {r['model']:58s} R2 {r['r2_mean']:.3f} \u00b1 {r['r2_sd']:.3f} "
                  f"({r['n_rows']} samples, {r['n_features']} predictors)")

    all_tiers = sorted(set(illumina.index) | set(full.index))
    x = np.arange(len(all_tiers))
    width = 0.35

    fig, ax = plt.subplots(figsize=(7.5, 5))
    err_kw = {"ecolor": "#333333", "capsize": 3, "elinewidth": 1}
    ax.bar(x - width / 2, [illumina.loc[t, "r2_mean"] if t in illumina.index else np.nan for t in all_tiers],
           width, yerr=[illumina.loc[t, "r2_sd"] if t in illumina.index else np.nan for t in all_tiers],
           label=f"Illumina-only ({illumina['folds'].iloc[0].count(';') + 1} studies)", color=CLR_GROUPED,
           edgecolor="black", linewidth=0.6, error_kw=err_kw)
    ax.bar(x + width / 2, [full.loc[t, "r2_mean"] if t in full.index else np.nan for t in all_tiers],
           width, yerr=[full.loc[t, "r2_sd"] if t in full.index else np.nan for t in all_tiers],
           label=f"Full matrix, cross-platform ({full['folds'].iloc[0].count(';') + 1} studies)",
           color=CLR_ACCENT, edgecolor="black", linewidth=0.6, error_kw=err_kw)
    for i, t in enumerate(all_tiers):
        if t not in full.index:
            ax.text(x[i] + width / 2, -0.08, "not run:\ncompounds\nmeasured in\nIllumina\nstudies only",
                    ha="center", va="top", fontsize=6.5, color="#555555", style="italic")

    ax.set_xlim(x[0] - 0.6, x[-1] + 0.6)  # keeps the slot of a tier that was not run inside the axes
    ax.axhline(0, color="#333333", linewidth=1.0)
    ax.set_ylabel("LOSO $R^2$ (mean \u00b1 SD across held-out studies)", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[t] for t in all_tiers])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=True, edgecolor="none")

    plt.title("Figure S1. Feature Curation Does Not Rescue Cross-Study Generalization (Table 4)",
              loc="left", fontweight="bold", pad=12, fontsize=9.5)
    plt.tight_layout()
    fi.save(fig, "Figure_S1_Substrate_Ladder")
    plt.close()


if __name__ == "__main__":
    generate_figure_s1()
