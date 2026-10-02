from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, finish


def draw(frame, kind):
    if kind == "heatmap":
        fig, ax = plt.subplots(figsize=(8, 5))
        image = ax.imshow(frame.to_numpy(float), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
        ax.set(xlabel="Tumor samples", ylabel="Signature genes", xticks=[])
        fig.colorbar(image, ax=ax, label="Expression z-score")
    elif kind == "survival":
        fig, ax = plt.subplots(figsize=(4, 3))
        for group, part in frame.groupby("group", sort=False):
            line, = ax.step(part.time, part.survival, where="post", label=group)
            censored = part.loc[part.censored > 0]
            ax.scatter(censored.time, censored.survival, marker="+", color=line.get_color())
        ax.set(xlabel="Time", ylabel="Progression-free survival", ylim=(0, 1.05))
        ax.legend(frameon=False)
    elif kind == "pdx":
        summary = frame.groupby("group").value.agg(["mean", "sem"]).reindex(["VS", "VR", "PS", "PR"])
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.bar(summary.index, summary["mean"], yerr=summary["sem"], capsize=3, color=["#aaaaaa", "#666666", "#dc9b26", "#b3423c"])
        ax.set_ylabel("Mean signature z-score ± SEM")
    else:
        genes = frame.gene.unique()
        fig, axes = plt.subplots(1, len(genes), figsize=(2 * len(genes), 3), squeeze=False)
        for ax, gene in zip(axes.flat, genes):
            part = frame.loc[frame.gene == gene]
            arrays = [part.loc[part.group == group, "value"].to_numpy(float) for group in ["Sensitive", "Resistant"]]
            ax.boxplot(arrays, tick_labels=["S", "R"], showfliers=False)
            for x, values in enumerate(arrays, 1):
                ax.scatter(x + np.linspace(-.12, .12, len(values)), values, s=9, color="black")
            ax.set(title=gene, ylabel="Normalized reads")
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--kind", choices=["heatmap", "pdx", "baseline", "survival"], required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    frame = pd.read_csv(a.input, sep="\t", index_col=0) if a.kind == "heatmap" else read_table(a.input)
    finish(draw(frame, a.kind), a.output)
