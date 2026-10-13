# figure_s4_duber_scatter.py
# Figure S4: actual vs. predicted caproate for the Duber et al. (2025) B1 -> B2
# transfer of Table 5 (Paper_ID Duber_2024). The model is refitted here exactly
# as 05 fits it (core set, training-reactor-only scaling, the shared settings
# with no subsampling), so the predictions are 05's. The script stops unless
# its R2, RMSE, Spearman rho and sample counts equal those in 05's run of
# record (05_table5_reactor_pairs.csv). Rho is on sampling-day means, as in
# Table 5.
# Run it with run_figures.py, which points it at the matrix and that file.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost as xgb

import figure_inputs as fi
from pipeline_settings import (CAT_COLS, CORE_CHEM, OPS, REACTOR_PAIRS, TARGET, XGB_PARAMS, load_matrix,
                               present, r2_rmse, sampling_units, scale, spearman_by_unit)

plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["mathtext.fontset"] = "custom"  # R2 and rho in the text font (Arial), not DejaVu Sans
plt.rcParams["mathtext.cal"] = "sans"  # unused; avoids a warning where no cursive font is installed
plt.rcParams["font.size"] = 9
plt.rcParams["axes.edgecolor"] = "#333333"
plt.rcParams["axes.linewidth"] = 0.8

CLR_NAIVE = "#E64B35"
PAPER = "Duber_2024"


def generate_figure_s4():
    """Figure S4 -- the real held-out predictions of Table 5's Duber 2025 pair."""
    df, genus_cols = load_matrix(fi.matrix_path())
    units = sampling_units(df)
    valid_features = sorted(set(present(OPS + CORE_CHEM + CAT_COLS, df) + genus_cols))  # as in 05
    br_train, br_test = REACTOR_PAIRS[PAPER]
    train = df[(df["Paper_ID"] == PAPER) & (df["BIOREACTOR"] == br_train)]
    test = df[(df["Paper_ID"] == PAPER) & (df["BIOREACTOR"] == br_test)]

    X_train, X_test = scale(train[valid_features], test[valid_features])
    model = xgb.XGBRegressor(**XGB_PARAMS).fit(X_train, train[TARGET])
    y_test, y_pred = test[TARGET], model.predict(X_test)
    r2, rmse = r2_rmse(y_test, y_pred)
    rho, n_units = spearman_by_unit(y_test, y_pred, units[test.index])

    t5 = fi.one_row(pd.read_csv(fi.results_file("05", "05_table5_reactor_pairs.csv")),
                    Paper=PAPER, Train_Reactor=br_train, Test_Reactor=br_test)
    for name, ours, theirs in [("R2", r2, t5["R2"]), ("RMSE", rmse, t5["RMSE (mM C)"]),
                               ("Spearman rho", rho, t5["Spearman_rho_unit_means"]),
                               ("train samples", len(train), t5["Train_n"]), ("test samples", len(test), t5["Test_n"]),
                               ("test sampling days", n_units, t5["Test_sampling_units"])]:
        if not np.isclose(ours, theirs, rtol=0, atol=1e-9):
            raise SystemExit(f"{name}: {ours} here, {theirs} in 05_table5_reactor_pairs.csv.")
    print(f"\n[Figure S4] {PAPER} {br_train} ({len(train)} samples) -> {br_test} ({len(test)} samples, "
          f"{n_units} sampling days): R2 = {r2:.3f}, RMSE = {rmse:.1f} mM C, Spearman rho (day means) = "
          f"{rho:.3f}; equal to Table 5 (05_table5_reactor_pairs.csv)")
    print(pd.DataFrame({"Sample_ID": test["Sample_ID"], "DAY": test["DAY"], "measured": y_test,
                        "predicted": y_pred}).to_string(index=False, float_format=lambda x: f"{x:.1f}"))

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(y_test, y_pred, color=CLR_NAIVE, alpha=0.85, edgecolor="black", s=45, linewidth=0.6)

    lims = [min(y_test.min(), y_pred.min()) - 10, max(y_test.max(), y_pred.max()) + 10]
    ax.plot(lims, lims, "k--", linewidth=1.0, label="Perfect prediction ($y = x$)")

    ax.set_xlabel("Actual Caproate (mM C)", fontweight="bold")
    ax.set_ylabel("Predicted Caproate (mM C)", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper left", frameon=True, edgecolor="none", fontsize=8)
    ax.text(0.97, 0.03, f"$R^2 = {r2:.3f}$ ({len(test)} samples)\n"
                        f"Spearman $\\rho = {rho:.3f}$ ({n_units} sampling-day means)", transform=ax.transAxes,
            fontsize=8.5, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#F8F9FA", edgecolor="#CCCCCC"))

    plt.title("Figure S4. Rank-Inverted Prediction (Duber et al. 2025, B1 \u2192 B2)",
              loc="left", fontweight="bold", pad=12, fontsize=9.5)
    plt.tight_layout()
    fi.save(fig, "Figure_S4_Duber_Scatter")
    plt.close()


if __name__ == "__main__":
    generate_figure_s4()
