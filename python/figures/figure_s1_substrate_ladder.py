import os
import numpy as np
import matplotlib.pyplot as plt

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
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

CLR_GROUPED = "#3C5488"
CLR_ACCENT = "#00A087"


def generate_figure_s1():
    """Figure S1 -- four-tier substrate ladder. Values are the same verified
    scalars already in manuscript Table 4 (CSTR-only and full-matrix LOSO R2
    per tier) -- a legitimate use of hardcoded values, same as main Figures 1
    and 3, since these are already-confirmed summary statistics, not
    per-sample data that could be simulated."""
    labels = [
        "Tier 1:\nOps + Genera",
        "Tier 2:\n+ Foundational Feeds",
        "Tier 3:\n+ Intermediates",
        "Tier 4:\n+ Downstream Proxies",
    ]
    cstr_r2 = [-0.163, -0.198, -0.376, -0.425]   # Table 4, CSTR-only column
    full_r2 = [-0.775, -0.923, -0.320, -0.095]   # Table 4, full-matrix column

    x = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.bar(x - width / 2, cstr_r2, width, label="CSTR-only universe (Illumina)", color=CLR_GROUPED, edgecolor="black", linewidth=0.6)
    ax.bar(x + width / 2, full_r2, width, label="Full matrix (cross-platform)", color=CLR_ACCENT, edgecolor="black", linewidth=0.6)

    ax.axhline(0, color="#333333", linewidth=1.0)
    ax.set_ylabel("LOSO $R^2$", fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower right", frameon=True, edgecolor="none")

    plt.title("Figure S1. Feature Curation Does Not Rescue Cross-Study Generalization (Table 4)",
              loc="left", fontweight="bold", pad=12, fontsize=9.5)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_S1_Substrate_Ladder.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_S1_Substrate_Ladder.pdf"), bbox_inches="tight")
    plt.close()
    print(f"[OK] Figure S1 saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    generate_figure_s1()
