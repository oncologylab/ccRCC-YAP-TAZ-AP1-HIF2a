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
    if kind == "clinical":
        table = frame.pivot(index="group", columns="category", values="percent").fillna(0)
        fig, ax = plt.subplots(figsize=(4, 3))
        bottom = np.zeros(len(table))
        for column in table:
            ax.bar(table.index, table[column], bottom=bottom, label=column)
            bottom += table[column]
        ax.set_ylabel("Patients (%)")
        ax.legend(frameon=False, bbox_to_anchor=(1.02, 1))
    elif kind == "heatmap":
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(frame.to_numpy(float), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
        ax.set(xticks=range(len(frame.columns)), xticklabels=frame.columns, ylabel="Clustered genes")
        ax.tick_params(axis="x", rotation=90)
        fig.colorbar(im, ax=ax, label="Expression z-score")
    else:
        clusters = frame.cluster.unique()
        fig, axes = plt.subplots(1, len(clusters), figsize=(4 * len(clusters), 4), squeeze=False)
        for ax, cluster in zip(axes.flat, clusters):
            part = frame.loc[frame.cluster == cluster].nsmallest(10, "pvalue").sort_values("pvalue", ascending=False)
            if not part.pvalue.between(0, 1, inclusive="right").all():
                raise ValueError("Pathway P values must be in (0, 1]")
            ax.barh(part.pathway, -np.log10(part.pvalue), color="#577aa7")
            ax.set(title=f"Cluster {cluster}", xlabel="-log10 P")
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--kind", choices=["clinical", "heatmap", "pathways"], required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    frame = pd.read_csv(a.input, sep="\t", index_col=0) if a.kind == "heatmap" else read_table(a.input)
    finish(draw(frame, a.kind), a.output)
