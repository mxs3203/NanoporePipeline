#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def jitter_points(n, width=0.25, seed=42):
    """Return n jitter offsets centered at 0 for a beeswarm-like scatter."""
    rng = np.random.default_rng(seed)
    return rng.uniform(-width, width, size=n)


def main():
    ap = argparse.ArgumentParser(
        description="Create bar + boxplot (with beeswarm) by biomarker_group from a methylation CSV."
    )
    ap.add_argument("csv", help="Input CSV with columns including biomarker_group, weight, and gene.")
    ap.add_argument("-o", "--out", default="biomarker_groups.png",
                    help="Output image filename (png/svg/pdf). Default: biomarker_groups.PDF")
    ap.add_argument("--only-biomarkers", action="store_true",
                    help="If set, keep only rows where is_biomarker == True.")
    ap.add_argument("--label-threshold", type=float, default=40.0,
                    help="Add gene labels to points with weight >= this threshold. Default: 40.0")
    args = ap.parse_args()

    # Load
    df = pd.read_csv(args.csv)

    # Basic column checks and normalization
    required_cols = {"biomarker_group", "weight", "gene"}
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Found columns: {list(df.columns)}")

    # Optional filter to biomarker rows only
    if args.only_biomarkers and "is_biomarker" in df.columns:
        df = df[df["is_biomarker"] == True].copy()

    # Coerce numerics
    df["weight"] = pd.to_numeric(df["weight"], errors="coerce")
    df["biomarker_group"] = df["biomarker_group"].astype(str).str.strip()
    df["gene"] = df["gene"].astype(str).str.strip()

    # Drop rows with missing group
    df = df[~df["biomarker_group"].isna() & (df["biomarker_group"] != "") & (df["biomarker_group"] != "nan")]
    df_box = df.dropna(subset=["weight"]).copy()

    # Order groups by count (descending)
    counts = df.groupby("biomarker_group", dropna=False).size().sort_values(ascending=False)
    groups = counts.index.tolist()

    # Prepare data for boxplot
    data_by_group = [df_box.loc[df_box["biomarker_group"] == g, "weight"].values for g in groups]

    # Create figure
    fig, (ax_bar, ax_box) = plt.subplots(2, 1, figsize=(12, 8), sharex=True,
                                         gridspec_kw={"height_ratios": [1, 1.4]})
    fig.suptitle("Biomarker Groups: Counts and Methylation Weight Distribution", y=0.98)

    # --- Top: Bar chart of counts ---
    x = np.arange(len(groups))
    ax_bar.bar(x, counts.values)
    ax_bar.set_ylabel("Count of genes")
    if len(counts):
        ax_bar.set_ylim(0, max(counts.values) * 1.15)
    for xi, yi in zip(x, counts.values):
        ax_bar.text(xi, yi, str(int(yi)), ha="center", va="bottom", fontsize=9)

    # --- Bottom: Boxplot + beeswarm ---
    # Boxplot (no outliers to keep scatter readable)
    ax_box.boxplot(
        data_by_group,
        positions=x,
        widths=0.6,
        patch_artist=False,
        showfliers=False,
        vert=True,
    )
    ax_box.set_ylabel("Methylation weight (%)")

    # Beeswarm + labels for points above threshold
    for xi, g in zip(x, groups):
        rows_g = df_box.loc[df_box["biomarker_group"] == g, ["weight", "gene"]]
        vals = rows_g["weight"].values
        genes = rows_g["gene"].values
        if len(vals) == 0:
            continue

        jitter = jitter_points(len(vals), width=0.22, seed=xi + 12345)  # per-group seed

        # Scatter the points
        ax_box.scatter(np.full_like(vals, xi, dtype=float) + jitter, vals, s=28, alpha=0.85)

        # Annotate points above threshold
        mask = vals >= args.label_threshold
        for j_off, v, gene_name in zip(jitter[mask], vals[mask], genes[mask]):
            # Slight vertical offset to keep labels legible
            ax_box.annotate(
                gene_name,
                xy=(xi + j_off, v),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=8,
                rotation=0,
                clip_on=True,
            )

    # X-axis labels
    ax_box.set_xticks(x)
    ax_box.set_xticklabels(groups, rotation=45, ha="right")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    print(f"Saved figure to: {out_path}")


if __name__ == "__main__":
    main()
