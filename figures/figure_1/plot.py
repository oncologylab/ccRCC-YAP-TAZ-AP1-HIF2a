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
    if kind == "cna":
        labels = frame.cohort + ": " + frame.gene
        table = frame.assign(label=labels).pivot(index="label", columns="status", values="percent").fillna(0)
        fig, ax = plt.subplots(figsize=(7, 3))
        bottom = np.zeros(len(table))
        for status, color in [("Loss", "#3b6fb6"), ("Neutral", "#cccccc"), ("Gain", "#c94b4b")]:
            values = table[status] if status in table else np.zeros(len(table))
            ax.bar(table.index, values, bottom=bottom, label=status, color=color)
            bottom += values
        ax.set_ylabel("Tumors (%)")
        ax.tick_params(axis="x", rotation=45)
        ax.legend(frameon=False)
    elif kind == "survival":
        fig, ax = plt.subplots(figsize=(4, 3))
        for group, part in frame.groupby("group", sort=False):
            ax.step(part.time, part.survival, where="post", label=group)
            censored = part.loc[part.censored > 0]
            ax.scatter(censored.time, censored.survival, marker="+", s=18)
        ax.set(xlabel="Time", ylabel="Survival probability", ylim=(0, 1.05))
        ax.legend(frameon=False)
    else:
        frame = frame.copy()
        frame["panel"] = frame.gene if "cohort" not in frame else frame.cohort + ": " + frame.gene
        panels = frame.panel.unique()
        fig, axes = plt.subplots(1, len(panels), figsize=(3 * len(panels), 3), squeeze=False)
        for ax, panel in zip(axes.flat, panels):
            part = frame.loc[frame.panel == panel]
            groups = [g for g in ["Neutral", "Loss"] if g in set(part.group)]
            for x, group in enumerate(groups):
                y = part.loc[part.group == group, "value"].to_numpy(float)
                ax.scatter(x + np.linspace(-.12, .12, len(y)), y, s=8, alpha=.45)
                ax.errorbar(x, y.mean(), y.std(ddof=1), color="black", capsize=4, fmt="_")
            ax.set(xticks=range(len(groups)), xticklabels=groups, title=panel, ylabel="Expression z-score")
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--kind", choices=["cna", "distribution", "survival"], required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    finish(draw(read_table(a.input), a.kind), a.output)
