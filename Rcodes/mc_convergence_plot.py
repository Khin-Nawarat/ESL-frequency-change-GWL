# Author: Khin Nawarat  (2026, @IHE Delft)
# Demonstrate whether 1000 realizations are sufficient for the 200-year event estimates (Figure 4)


from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


BASE_DIR = Path(r'D:\Khin Nawarat\ESLs at different GWLs\0_SD_runs') # change to your local directory
OUTPUT_DIR = BASE_DIR / "Rdatasets"
CSV_PATH = OUTPUT_DIR / "mc_convergence_summary.csv"
FIG_PATH = BASE_DIR / "pics" / "mc_convergence_plot.png"


if not CSV_PATH.exists():
    raise FileNotFoundError(f"Cannot find summary CSV at {CSV_PATH}")

# Read the convergence summary and use 1000 / 5000 sample sizes only.
df = pd.read_csv(CSV_PATH)
df = df[df["nsample"].isin([1000, 5000])]

# Calculate percentage differences for all three percentiles
# All as percentage of the present-day 50th percentile (p50_ref)
df["percentage_p05"] = (df["abs_diff_p05"] / df["p50_ref"]) * 100
df["percentage_p50"] = (df["abs_diff_p50"] / df["p50_ref"]) * 100
df["percentage_p95"] = (df["abs_diff_p95"] / df["p50_ref"]) * 100

df = df[["dataset", "SLRID", "nsample", "percentage_p05", "percentage_p50", "percentage_p95"]]

order = [
    "Rasmussen_1000", "Rasmussen_5000",
    "Vousdoukas_1000", "Vousdoukas_5000",
    "Kirezci_1000", "Kirezci_5000"
]

# Keep dataset order consistent
df["dataset"] = pd.Categorical(df["dataset"], categories=["Rasmussen", "Vousdoukas", "Kirezci"], ordered=True)
df["x_group"] = pd.Categorical(
    df["dataset"].astype(str) + "_" + df["nsample"].astype(str),
    categories=order,
    ordered=True
)

# Reshape to long format for easier grouping by percentile
df_long = df.melt(
    id_vars=["dataset", "SLRID", "nsample", "x_group"],
    value_vars=["percentage_p05", "percentage_p50", "percentage_p95"],
    var_name="percentile",
    value_name="percentage"
)

percentile_map = {
    "percentage_p05": "5th percentile",
    "percentage_p50": "50th percentile",
    "percentage_p95": "95th percentile"
}
df_long["percentile"] = df_long["percentile"].map(percentile_map)

summary = df_long.groupby(["percentile", "dataset", "nsample", "x_group"], observed=True)["percentage"].agg(
    q05=lambda x: x.quantile(0.05),
    q25=lambda x: x.quantile(0.25),
    q50=lambda x: x.quantile(0.50),
    q75=lambda x: x.quantile(0.75),
    q95=lambda x: x.quantile(0.95)
).reset_index()

def draw_plot(fig_path, y_lim):
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    for ax, percentile in zip(axes, ["5th percentile", "50th percentile", "95th percentile"]):
        subset = summary[summary["percentile"] == percentile]
        x = list(range(len(order)))
        for idx, row in subset.iterrows():
            xpos = order.index(row["x_group"])
            ax.vlines(xpos, row["q05"], row["q95"], color="black", linewidth=0.8)
            rect = Rectangle(
                (xpos - 0.18, row["q25"]),
                width=0.36,
                height=max(row["q75"] - row["q25"], 1e-9),
                facecolor="lightgray",
                edgecolor="black",
                linewidth=0.8,
                zorder=2
            )
            ax.add_patch(rect)
            ax.hlines(row["q50"], xpos - 0.18, xpos + 0.18, color="Black", linewidth=1.5, zorder=3)
            ax.hlines(row["q05"], xpos - 0.08, xpos + 0.08, color="black", linewidth=0.8)
            ax.hlines(row["q95"], xpos - 0.08, xpos + 0.08, color="black", linewidth=0.8)

        ax.set_title(percentile, fontsize=18, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels([label.split("_")[1] for label in order], fontsize=14, fontweight="bold")
        ax.set_xlim(-0.6, len(order) - 0.4)
        ax.set_ylim(*y_lim)
        ax.set_xlabel("Sample size", fontsize=18, fontweight="bold")
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(0.8)
        ax.spines["bottom"].set_linewidth(0.8)
        ax.tick_params(axis="x", length=0)
        ax.tick_params(axis="y", labelsize=14)
        for center, label in zip([0.5, 2.5, 4.5], ["Rasmussen", "Vousdoukas", "Kirezci"]):
            ax.text(
                center,
                0.9,
                label,
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=16,
                fontweight="bold"
            )

    axes[0].set_ylabel("Percentage of present-day\n median 200-year ESL", fontsize=18, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(fig_path, dpi=150)
    print(f"Saved plot to {fig_path}")

# Create the figure
draw_plot(FIG_PATH, (-0.5, 10))
