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
    if kind == "signal":
        panels = frame.comparison.unique()
        fig, axes = plt.subplots(1, len(panels), figsize=(3 * len(panels), 3), squeeze=False)
        for ax, panel in zip(axes.flat, panels):
            part = frame.loc[frame.comparison == panel]
            means = [part.control.mean(), part.treated.mean()]
            sd = [part.control.std(), part.treated.std()]
            ax.bar(["Ctrl", "KD"], means, yerr=sd, capsize=3, color=[".6", "#b24d4d"])
            ax.set(title=str(panel), ylabel="Mean H3K27ac signal")
    elif kind == "transcription-factors":
        table = frame.pivot(index="subset", columns="status", values="percent")
        fig, ax = plt.subplots(figsize=(4, 3))
        bottom = np.zeros(len(table))
        for status, color in [("Down", "#4a78ac"), ("NC", ".8"), ("Up", "#b45454")]:
            ax.bar(table.index, table[status], bottom=bottom, label=status, color=color)
            bottom += table[status]
        ax.set_ylabel("Transcription factors (%)")
        ax.legend(frameon=False)
    else:
        panels = frame.panel.unique()
        fig, axes = plt.subplots(1, len(panels), figsize=(3 * len(panels), 3), squeeze=False)
        for ax, panel in zip(axes.flat, panels):
            for group, part in frame.loc[frame.panel == panel].groupby("group", sort=False):
                part = part.sort_values("x")
                ax.plot(part.x, part.value, label=group)
            ax.set(title=str(panel), xlabel="bp from peak centre", ylabel="Signal")
            ax.legend(frameon=False, fontsize=7)
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--kind", choices=["signal", "transcription-factors", "profile"], required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    finish(draw(read_table(a.input), a.kind), a.output)
