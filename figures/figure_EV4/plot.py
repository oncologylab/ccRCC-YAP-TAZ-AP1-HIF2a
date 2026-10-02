from pathlib import Path
import argparse
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, finish


def draw(frame, kind):
    if kind == "occupancy":
        frame = frame.sort_values("n_regions")
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.barh(frame["class"], frame.n_regions, color="#78689a")
        ax.set_xlabel("Regions")
    elif kind == "heatmap":
        fig, ax = plt.subplots(figsize=(6, 5))
        image = ax.imshow(frame.to_numpy(float), aspect="auto", cmap="magma", interpolation="none")
        ax.set(xticks=range(len(frame.columns)), xticklabels=frame.columns, ylabel="Stratified peaks")
        ax.tick_params(axis="x", rotation=90)
        fig.colorbar(image, ax=ax, label="CUT&Tag signal")
    else:
        panels = frame.panel.unique()
        fig, axes = plt.subplots(1, len(panels), figsize=(3.2 * len(panels), 3), squeeze=False)
        for ax, panel in zip(axes.flat, panels):
            for group, part in frame.loc[frame.panel == panel].groupby("group", sort=False):
                part = part.sort_values("x")
                ax.plot(part.x, part.value, label=group)
            ax.set(title=str(panel), xlabel="Distance from centre (bp)", ylabel="Mean signal")
            ax.legend(frameon=False, fontsize=7)
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--kind", choices=["occupancy", "heatmap", "profiles"], required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    frame = pd.read_csv(a.input, sep="\t", index_col=0) if a.kind == "heatmap" else read_table(a.input)
    finish(draw(frame, a.kind), a.output)
