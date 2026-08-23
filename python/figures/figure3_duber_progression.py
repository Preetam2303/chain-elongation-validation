# figure3_duber_progression.py
# The three-step leakage-tightening illustration: the identical Duber et al.
# (2025) B1(n=30) -> B2(n=21) split evaluated three times as leakage controls
# were progressively tightened, with no change to the underlying samples.
# Values are the same verified scalars reported in Table 5 and Section 3.4
# (see python/pipeline/07a, 07c, and 05 for the scripts behind each of the
# three numbers respectively).

import os
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
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8


def generate_figure_3():
    stages = [
        "1. Unscaled Random Forest\n(+ Downstream Proxies)",
        "2. XGBoost\n(No Downstream Exclusion,\nUnconfirmed Scaling)",
        "3. Fully Leakage-Controlled\nXGBoost (Train-Only Scaling)",
    ]
    r2_values = [0.217, -0.228, -5.869]
    colors = ["#4DBBD5", "#E64B35", "#7E6148"]

    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    bars = ax.bar(stages, r2_values, color=colors, width=0.55, edgecolor="black", linewidth=0.8)

    ax.axhline(0, color="#333333", linestyle="-", linewidth=0.8)
    ax.set_ylabel(r"Transferability ($R^2$, Reactor B1 $\rightarrow$ B2)", fontweight="bold")
    ax.set_ylim(-7.0, 1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for bar in bars:
        yval = bar.get_height()
        offset = 0.2 if yval >= 0 else -0.4
        ax.text(bar.get_x() + bar.get_width() / 2, yval + offset, f"$R^2 = {yval:.3f}$",
                ha="center", va="bottom" if yval >= 0 else "top", fontweight="bold", fontsize=8.5)

    plt.title("Figure 3. Duber et al. (2025) Intra-Study Transfer Progression",
              loc="left", fontweight="bold", pad=12)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_3_Duber_Progression.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_3_Duber_Progression.pdf"), bbox_inches="tight")
    plt.close()
    print(f"[OK] Figure 3 saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    generate_figure_3()
