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
plt.rcParams["xtick.labelsize"] = 8
plt.rcParams["ytick.labelsize"] = 8.5
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

CLR_NAIVE = "#E64B35"
CLR_GROUPED = "#3C5488"

def generate_figure_s3():
    studies = [
        "Brodowski 2022 [Ext.]\n(B1 \u2192 B2)",
        "Brodowski 2025\n(B1 \u2192 B2)",
        "Duber 2022\n(R1 \u2192 R2)",
        "Prusak 2025\n(B1 \u2192 B2)",
        "Duber 2025\n(B1 \u2192 B2)",
    ]
    r2_vals =  [0.442, 0.400, 0.008, -0.140, -5.869]
    rho_vals = [0.324, 0.964, 1.000, 0.792, -0.569]

    x = np.arange(len(studies))
    width = 0.35

    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    bars_r2 = ax.bar(x - width / 2, r2_vals, width, label="$R^2$",
                     color=[CLR_GROUPED if v >= 0 else CLR_NAIVE for v in r2_vals],
                     edgecolor="black", linewidth=0.6)
    bars_rho = ax.bar(x + width / 2, rho_vals, width, label="Spearman $\\rho$", color="#8491B4",
                      edgecolor="black", linewidth=0.6, hatch="//")

    ax.axhline(0, color="#333333", linewidth=1.0)
    ax.set_ylabel("Value (transfer $R^2$ or rank correlation)", fontweight="bold", labelpad=10)
    ax.set_ylim(-7.0, 1.5)
    ax.set_xticks(x)
    ax.set_xticklabels(studies)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="lower left", frameon=True, edgecolor="none")

    for bar in bars_r2:
        yval = bar.get_height()
        offset = 0.2 if yval >= 0 else -0.35
        ax.text(bar.get_x() + bar.get_width() / 2, yval + offset, f"{yval:.3f}",
                ha="center", va="bottom" if yval >= 0 else "top", fontweight="bold", fontsize=7.5)
    for bar in bars_rho:
        yval = bar.get_height()
        offset = 0.2 if yval >= 0 else -0.35
        ax.text(bar.get_x() + bar.get_width() / 2, yval + offset, f"{yval:.3f}",
                ha="center", va="bottom" if yval >= 0 else "top", fontweight="bold", fontsize=7.5)

    plt.title("Figure S3 (optional). Intra-Study Reactor-Pair Transferability: $R^2$ vs. Rank Correlation (Table 5)",
              loc="left", fontweight="bold", pad=12, fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_S3_Intrastudy_Transfer.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure_S3_Intrastudy_Transfer.pdf"), bbox_inches="tight")
    plt.close()
    print(f"[OK] Figure S3 saved to {OUTPUT_DIR}")

if __name__ == "__main__":
    generate_figure_s3()
