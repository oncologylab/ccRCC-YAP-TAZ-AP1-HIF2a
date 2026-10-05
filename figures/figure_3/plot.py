from pathlib import Path
import argparse
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, finish


def draw(frame, kind, ylabel):
    if kind == "matched":
        panels = frame.panel.unique()
        fig, axes = plt.subplots(1, len(panels), figsize=(3.2 * len(panels), 3), squeeze=False)
        for ax, panel in zip(axes.flat, panels):
            part = frame.loc[frame.panel == panel]
            arrays = [part.loc[part.peak_status == group, "log2FC"] for group in ["bound", "unbound"]]
            ax.boxplot(arrays, tick_labels=["Co-occupied", "Background"], showfliers=False, whis=1.5)
            title = part.display_title.iloc[0] if "display_title" in part else str(panel)
            ax.set(title=title, ylabel="log2 fold change")
            ax.axhline(0, color="0.6", lw=.7)
    elif kind == "growth":
        fig, ax = plt.subplots(figsize=(4, 3))
        for group, part in frame.groupby("group", sort=False):
            part = part.sort_values("time")
            ax.errorbar(part.time, part["mean"], part["sem"], label=group, capsize=2, marker="o", markersize=3)
        ax.set(xlabel="Time", ylabel=ylabel)
        ax.legend(frameon=False, fontsize=7)
    elif kind == "overlap":
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.bar(frame.peak_set, frame.overlapping, label="Overlapping")
        ax.bar(frame.peak_set, frame.nonoverlapping, bottom=frame.overlapping, label="Non-overlapping", color=".7")
        ax.set_ylabel("Peaks")
        ax.legend(frameon=False)
    else:
        panels = frame.panel.unique()
        fig, axes = plt.subplots(1, len(panels), figsize=(3.3 * len(panels), 3), squeeze=False)
        for ax, panel in zip(axes.flat, panels):
            for group, part in frame.loc[frame.panel == panel].groupby("group", sort=False):
                part = part.sort_values("x")
                ax.plot(part.x, part.value, label=group)
            ax.set(title=str(panel), xlabel="bp from peak centre", ylabel=ylabel)
            ax.legend(frameon=False, fontsize=7)
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--kind", choices=["matched", "growth", "overlap", "profile"], required=True)
    p.add_argument("--ylabel", default="Signal")
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    finish(draw(read_table(a.input), a.kind, a.ylabel), a.output)
