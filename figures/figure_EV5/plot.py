from pathlib import Path
import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, finish


def draw(frame, kind):
    if kind == "composition":
        table = frame.pivot(index="stratum", columns="motif_pattern", values="percent").fillna(0)
        fig, ax = plt.subplots(figsize=(6, 3))
        bottom = np.zeros(len(table))
        for pattern in table:
            ax.bar(table.index, table[pattern], bottom=bottom, label=pattern)
            bottom += table[pattern]
        ax.set_ylabel("Peaks (%)")
        ax.legend(frameon=False, fontsize=7, bbox_to_anchor=(1.02, 1))
    elif kind == "distance":
        fig, ax = plt.subplots(figsize=(5, 3))
        finite = frame.loc[np.isfinite(frame.log10_distance_plus1)]
        bins = np.linspace(0, max(1., finite.log10_distance_plus1.max()), 33)
        for family, part in finite.groupby("family"):
            ax.hist(part.log10_distance_plus1, bins=bins, histtype="step", density=True, label=family, linewidth=1.5)
        ax.set(xlabel="log10(nearest midpoint distance + 1 bp)", ylabel="Density")
        ax.legend(frameon=False)
    elif kind == "volcano":
        fig, ax = plt.subplots(figsize=(5, 4))
        for status, color in [("NS", ".75"), ("Ctrl", "#476bb2"), ("KD", "#b84d4d")]:
            part = frame.loc[frame.status == status]
            ax.scatter(part.differential_score, part.minus_log10_padj, c=color, label=status, s=12)
        for row in frame.loc[frame.motif.str.contains("JUN|TEAD|HIF|ARNT", case=False, regex=True)].itertuples():
            ax.annotate(row.motif, (row.differential_score, row.minus_log10_padj), fontsize=6)
        ax.set(xlabel="Differential footprint score (KD − Ctrl)", ylabel="-log10 adjusted P")
        ax.legend(frameon=False)
    else:
        panels = frame.panel.unique()
        fig, axes = plt.subplots(1, len(panels), figsize=(3 * len(panels), 3), squeeze=False)
        for ax, panel in zip(axes.flat, panels):
            part = frame.loc[frame.panel == panel]
            groups = part.group.unique()
            for x, group in enumerate(groups):
                y = part.loc[part.group == group, "value"].to_numpy(float)
                ax.scatter(x + np.linspace(-.15, .15, len(y)), y, s=6, alpha=.4)
                ax.errorbar(x, y.mean(), y.std(ddof=1), fmt="_", color="black", capsize=3)
            ax.set(title=str(panel), xticks=range(len(groups)), xticklabels=groups, ylabel="Nuclear fluorescence intensity")
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--kind", choices=["composition", "distance", "volcano", "imaging"], required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    finish(draw(read_table(a.input), a.kind), a.output)
