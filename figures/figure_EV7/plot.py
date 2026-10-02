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
    if kind == "intersections":
        frame = frame.sort_values("n_genes")
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.barh(frame.membership, frame.n_genes, color="#926698")
        ax.set_xlabel("Downregulated genes")
    elif kind == "chromatin":
        strata = frame.stratum.unique()
        fig, axes = plt.subplots(1, len(strata), figsize=(3 * len(strata), 3), squeeze=False)
        for ax, stratum in zip(axes.flat, strata):
            part = frame.loc[frame.stratum == stratum]
            means = [part.control.mean(), part.shJUNs.mean()]
            sd = [part.control.std(), part.shJUNs.std()]
            ax.bar(["Ctrl", "shJUNs"], means, yerr=sd, capsize=3, color=[".6", "#916ca7"])
            ax.set(title=str(stratum), ylabel="Mean H3K27ac signal")
    else:
        subsets = frame.membership.unique()
        fig, axes = plt.subplots(1, len(subsets), figsize=(4.5 * len(subsets), 4), squeeze=False)
        for ax, subset in zip(axes.flat, subsets):
            part = frame.loc[frame.membership == subset].nsmallest(10, "pvalue").sort_values("pvalue", ascending=False)
            if not part.pvalue.between(0, 1, inclusive="right").all():
                raise ValueError("Pathway P values must be in (0, 1]")
            ax.barh(part.pathway, -np.log10(part.pvalue), color="#926698")
            ax.set(title=str(subset), xlabel="-log10 P")
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--kind", choices=["intersections", "chromatin", "pathways"], required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    finish(draw(read_table(a.input), a.kind), a.output)
