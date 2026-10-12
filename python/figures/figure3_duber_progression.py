# figure3_duber_progression.py
# The three-step leakage-tightening illustration: the identical Duber et al.
# (2025) B1 -> B2 split (Paper_ID Duber_2024) evaluated three times as leakage
# controls were progressively tightened, with no change to the underlying
# samples. Each value is read from the corrected run's run of record:
#   stage 1  07d_figure3_stage1.csv (unscaled random forest, broad predictor list)
#   stage 2  07c_figure3_stage2.csv (XGBoost, broad predictor list, no scaling)
#   stage 3  05_table5_reactor_pairs.csv (Table 5)
# The script stops if the three stages do not use the same train and test
# samples. Run it with run_figures.py, which points it at those files.

import matplotlib.pyplot as plt
import pandas as pd

import figure_inputs as fi

plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = 9
plt.rcParams["axes.labelsize"] = 10
plt.rcParams["axes.titlesize"] = 10
plt.rcParams["xtick.labelsize"] = 8.5
plt.rcParams["ytick.labelsize"] = 8.5
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

PAIR = dict(Paper="Duber_2024", Train_Reactor="B1", Test_Reactor="B2")


def generate_figure_3():
    stage1 = fi.one_row(pd.read_csv(fi.results_file("07d", "07d_figure3_stage1.csv")), **PAIR)
    stage2 = fi.one_row(pd.read_csv(fi.results_file("07c", "07c_figure3_stage2.csv")), **PAIR)
    stage3 = fi.one_row(pd.read_csv(fi.results_file("05", "05_table5_reactor_pairs.csv")), **PAIR)
    sizes = {(int(stage1["Train_Size"]), int(stage1["Test_Size"])),
             (int(stage2["Train_Size"]), int(stage2["Test_Size"])),
             (int(stage3["Train_n"]), int(stage3["Test_n"]))}
    if len(sizes) != 1:
        raise SystemExit(f"The three stages do not use the same split: {sizes}")
    n_train, n_test = sizes.pop()

    stages = [
        "1. Unscaled Random Forest\n(+ Downstream Proxies)",
        "2. XGBoost\n(No Downstream Exclusion,\nUnconfirmed Scaling)",
        "3. Fully Leakage-Controlled\nXGBoost (Train-Only Scaling)",
    ]
    r2_values = [float(stage1["R2"]), float(stage2["R2"]), float(stage3["R2"])]
    colors = ["#4DBBD5", "#E64B35", "#7E6148"]
    print(f"\n[Figure 3] Duber_2024 B1 ({n_train} samples) -> B2 ({n_test} samples):")
    for s, v in zip(["Stage 1 (07d)", "Stage 2 (07c)", "Stage 3 (05, Table 5)"], r2_values):
        print(f"  {s}: R2 = {v:.3f} ({v!r})")

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    bars = ax.bar(stages, r2_values, color=colors, width=0.55, edgecolor="black", linewidth=0.8)

    ax.axhline(0, color="#333333", linestyle="-", linewidth=0.8)
    ax.set_ylabel(r"Transferability ($R^2$, Reactor B1 $\rightarrow$ B2)", fontweight="bold")
    ax.set_ylim(min(r2_values) - 1.1, max(1.5, max(r2_values) + 0.8))
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
    fi.save(fig, "Figure_3_Duber_Progression")
    plt.close()


if __name__ == "__main__":
    generate_figure_3()
