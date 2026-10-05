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
    panels = frame.panel.unique()
    fig, axes = plt.subplots(1, len(panels), figsize=(4.5 * len(panels), 3), squeeze=False)
    for ax, panel in zip(axes.flat, panels):
        part = frame.loc[frame.panel == panel]
        antibodies = part.antibody.unique()
        conditions = part.condition.unique()
        width = .8 / len(conditions)
        for j, condition in enumerate(conditions):
            for i, antibody in enumerate(antibodies):
                y = part.loc[(part.condition == condition) & (part.antibody == antibody), "percent_input"].to_numpy(float)
                if not len(y):
                    continue
                x = i + (j - (len(conditions) - 1) / 2) * width
                ax.bar(x, y.mean(), width=width, label=condition if i == 0 else None, color=plt.get_cmap("tab10")(j))
                if len(y) > 1:
                    ax.errorbar(x, y.mean(), y.std(ddof=1) / np.sqrt(len(y)), color="black", capsize=2)
                ax.scatter(x + np.linspace(-width / 6, width / 6, len(y)), y, s=10, color="black")
        ax.set(title=str(panel), xticks=range(len(antibodies)), xticklabels=antibodies, ylabel="Input enrichment (%)")
        ax.legend(frameon=False, fontsize=7)
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    finish(draw(read_table(a.input, ["panel", "antibody", "condition", "replicate", "percent_input"])), a.output)
