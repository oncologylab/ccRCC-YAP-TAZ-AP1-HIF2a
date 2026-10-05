from pathlib import Path
import argparse
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, finish


def draw(frame, locus, features=None):
    frame = frame.loc[frame.locus == locus]
    if frame.empty:
        raise ValueError(f"Missing locus: {locus}")
    tracks = frame.track.unique()
    nrows = len(tracks) + (features is not None)
    fig, axes = plt.subplots(nrows, 1, figsize=(9, max(3, nrows * .85)), sharex=True, squeeze=False)
    for ax, track in zip(axes.flat, tracks):
        part = frame.loc[frame.track == track].sort_values("x")
        ax.fill_between(part.x, part.value, color="#526f96", linewidth=0)
        ax.plot(part.x, part.value, color="#31486c", lw=.45)
        scale = frame.loc[frame.scale_group == part.scale_group.iloc[0], "value"].dropna()
        low, high = min(0, scale.min()), max(0, scale.max())
        ax.set_ylim(low, high if high > low else low + 1)
        ax.set_ylabel(track, rotation=0, ha="right", va="center", fontsize=7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="y", labelsize=6)
    if features is not None:
        ax = axes[-1, 0]
        selected = features.loc[features.locus == locus]
        for row in selected.itertuples(index=False):
            ax.axvspan(row.start, row.end, ymin=.25, ymax=.65, color="#8f6374", alpha=.6)
            ax.text((row.start + row.end) / 2, .75, row.label, ha="center", fontsize=7)
        ax.set(ylim=(0, 1), yticks=[])
        ax.spines[["top", "right", "left"]].set_visible(False)
    axes[0, 0].set_title(locus)
    axes[-1, 0].set_xlabel(f"{frame.chrom.iloc[0]} position (Mb)")
    axes[-1, 0].xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x / 1e6:.3f}"))
    return fig


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--locus", required=True)
    p.add_argument("--features", type=Path)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    features = read_table(a.features, ["locus", "start", "end", "label"]) if a.features else None
    finish(draw(read_table(a.input, ["locus", "track", "scale_group", "chrom", "x", "value"]), a.locus, features), a.output)
