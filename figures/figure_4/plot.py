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
    if kind == "counts":
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.bar(frame["class"], frame.n_regions, color="#6d83b3")
        ax.set_ylabel("Regions")
        ax.tick_params(axis="x", rotation=35)
        return fig
    panels = frame.panel.unique()
    fig, axes = plt.subplots(1, len(panels), figsize=(3.4 * len(panels), 3), squeeze=False)
    for ax, panel in zip(axes.flat, panels):
        part = frame.loc[frame.panel == panel]
        if kind == "imaging":
            groups = part.group.unique()
            for x, group in enumerate(groups):
                y = part.loc[part.group == group, "value"].to_numpy(float)
                ax.scatter(x + np.linspace(-.15, .15, len(y)), y, s=7, alpha=.5)
                ax.errorbar(x, y.mean(), y.std(ddof=1), fmt="_", color="black", capsize=4)
            ax.set(xticks=range(len(groups)), xticklabels=groups, ylabel="Overlapping nuclear area (%)")
        else:
            for condition, values in part.groupby("condition", sort=False):
                values = values.sort_values("position")
                ax.plot(values.position, values["mean"], label=condition)
            ax.axvline(0, color=".7", lw=.5)
            ax.set(xlabel="bp from motif centre", ylabel="Corrected Tn5 signal")
            ax.legend(frameon=False, fontsize=7)
        ax.set_title(str(panel))
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--kind", choices=["imaging", "footprints", "counts"], required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    finish(draw(read_table(a.input), a.kind), a.output)
