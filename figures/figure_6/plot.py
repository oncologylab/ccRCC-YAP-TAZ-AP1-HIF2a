from pathlib import Path
import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, finish


def draw(frame, panel):
    if panel == "C":
        clusters = sorted(frame.cluster.unique())
        fig, axes = plt.subplots(1, len(clusters), figsize=(2.4 * len(clusters), 3), squeeze=False)
        for ax, cluster in zip(axes.flat, clusters):
            part = frame.loc[frame.cluster == cluster]
            ax.boxplot([part.loc[part.group == group, "value"] for group in ["Ctrl", "shJUN"]],
                       tick_labels=["Ctrl", "shJUNs"], showfliers=False, whis=1.5)
            ax.set(title=f"C{cluster}", ylabel="Mean per-gene z-score")
    elif panel == "D":
        genes = frame.gene.unique()
        fig, axes = plt.subplots(1, len(genes), figsize=(2.7 * len(genes), 3), squeeze=False)
        for ax, gene in zip(axes.flat, genes):
            part = frame.loc[frame.gene == gene]
            for x, contrast in enumerate(["shYT", "shHIF2a", "shJUNs"]):
                values = part.loc[part.contrast == contrast, "displayed_log2FC"].to_numpy(float)
                ax.bar(x, values.mean(), yerr=values.std(ddof=1), color=["#d5a52b", "#b14b4b", "#667faf"][x], capsize=2)
                ax.scatter(x + np.linspace(-.08, .08, len(values)), values, s=13, color="black")
            ax.set(title=gene, xticks=range(3), xticklabels=["shYT", "shHIF2a", "shJUNs"], ylabel="log2 fold change")
            ax.tick_params(axis="x", rotation=35)
    else:
        fig, ax = plt.subplots(figsize=(3, 3))
        for x, group in enumerate(["Ctrl", "shJUNs"]):
            values = frame.loc[frame.group == group, "value"].to_numpy(float)
            ax.scatter(x + np.linspace(-.14, .14, len(values)), values, s=7, alpha=.5)
            ax.errorbar(x, values.mean(), values.std(ddof=1), fmt="_", color="black", capsize=4)
        ax.set(xticks=[0, 1], xticklabels=["Ctrl", "shJUNs"], ylabel="Overlapping nuclear area (%)")
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--panel", choices=["B", "C", "D"], required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    finish(draw(read_table(a.input), a.panel), a.output)
