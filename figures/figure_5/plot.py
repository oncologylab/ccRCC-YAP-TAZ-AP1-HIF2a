from pathlib import Path
import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, finish


def draw(frame):
    if not np.isfinite(frame.relative_binding).all():
        raise ValueError("Non-finite binding measurement")
    if frame.duplicated(["panel", "target", "condition", "replicate"]).any():
        raise ValueError("Duplicate measurement")
    panels = frame[["panel", "target"]].drop_duplicates()
    fig, axes = plt.subplots(1, len(panels), figsize=(3.6 * len(panels), 3), squeeze=False)
    for ax, (panel, target) in zip(axes.flat, panels.itertuples(index=False, name=None)):
        part = frame.loc[(frame.panel == panel) & (frame.target == target)]
        groups = part.condition.unique()
        for x, group in enumerate(groups):
            y = part.loc[part.condition == group, "relative_binding"].to_numpy(float)
            ax.bar(x, y.mean(), color=".8", width=.6)
            if len(y) > 1:
                ax.errorbar(x, y.mean(), y.std(ddof=1), capsize=3, color="black")
            ax.scatter(x + np.linspace(-.1, .1, len(y)), y, color="black", s=15)
        ax.set(title=f"{panel}: {target}", xticks=range(len(groups)), xticklabels=groups, ylabel="Relative IP/input")
        ax.tick_params(axis="x", rotation=35)
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    finish(draw(read_table(a.input, ["panel", "target", "condition", "replicate", "relative_binding"])), a.output)
