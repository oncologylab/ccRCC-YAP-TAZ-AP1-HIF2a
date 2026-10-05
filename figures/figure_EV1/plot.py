from pathlib import Path
import argparse
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import linregress

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, finish


def draw(frame, kind):
    fig, ax = plt.subplots(figsize=(5, 4))
    if kind == "methylation":
        ax.scatter(frame.methylation_beta, frame.expression_z, s=10, alpha=.5, color="#41658a")
        fit = linregress(frame.methylation_beta, frame.expression_z)
        x = frame.methylation_beta.sort_values()
        ax.plot(x, fit.intercept + fit.slope * x, color="black", lw=.8)
        ax.set(xlabel="SAV1 methylation beta", ylabel="SAV1 mRNA z-score")
    elif kind == "signature":
        frame = frame.sort_values("rank")
        ax.barh(frame.gene, frame.log2FC, color="#aa5555")
        ax.invert_yaxis()
        ax.set_xlabel("log2 fold change")
        fig.set_size_inches(5, max(3, len(frame) * .18))
    else:
        im = ax.imshow(frame.to_numpy(float), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
        ax.set(xlabel="Tumors", ylabel="Signature genes", xticks=[])
        fig.colorbar(im, ax=ax, label="Expression z-score")
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--kind", choices=["methylation", "signature", "heatmap"], required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    frame = pd.read_csv(a.input, sep="\t", index_col=0) if a.kind == "heatmap" else read_table(a.input)
    finish(draw(frame, a.kind), a.output)
