from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def read_table(path, required=()):
    path = Path(path)
    frame = pd.read_csv(path, sep="," if path.suffix == ".csv" else "\t")
    missing = set(required) - set(frame.columns)
    if missing or frame.empty:
        raise ValueError(f"Empty input or missing columns: {sorted(missing)}")
    return frame


def new_output(path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=False)
    return path


def save_table(frame, path):
    frame.to_csv(path, sep="\t", index=False)


def finish(fig, path):
    import matplotlib.pyplot as plt
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def bh(values):
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Non-finite P values")
    return false_discovery_control(values, method="bh")


def numeric(frame, columns):
    values = frame[list(columns)].to_numpy(float)
    if not np.isfinite(values).all():
        raise ValueError("Non-finite measurements")
    return values


def read_bed(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith(("#", "track", "browser")):
            continue
        fields = line.split("\t")
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        if start < 0 or end <= start:
            raise ValueError("Invalid BED interval")
        rows.append((chrom, start, end))
    if not rows:
        raise ValueError("Empty BED file")
    return pd.DataFrame(rows, columns=["chrom", "start", "end"]).drop_duplicates().reset_index(drop=True)


def overlap_flags(regions, sites):
    flags = np.zeros(len(regions), dtype=bool)
    regions = regions.reset_index(drop=True)
    parts = {chrom: part.sort_values("start") for chrom, part in sites.groupby("chrom")}
    for chrom, part in regions.groupby("chrom"):
        target = parts.get(chrom)
        if target is None:
            continue
        starts = target.start.to_numpy()
        ends = np.maximum.accumulate(target.end.to_numpy())
        stops = np.searchsorted(starts, part.end.to_numpy(), side="left")
        flags[part.index] = (stops > 0) & (ends[np.maximum(stops - 1, 0)] > part.start.to_numpy())
    return flags
