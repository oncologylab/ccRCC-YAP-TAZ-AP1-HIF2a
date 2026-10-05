from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, new_output, save_table


def locus_profiles(tracks, loci, bins, missing):
    import pyBigWig
    if bins < 1 or tracks.track.duplicated().any() or loci.locus.duplicated().any():
        raise ValueError("Positive bins and unique track/locus names required")
    records = []
    for track in tracks.itertuples(index=False):
        with pyBigWig.open(str(track.path)) as signal:
            for locus in loci.itertuples(index=False):
                start, end = int(locus.start), int(locus.end)
                length = signal.chroms().get(locus.chrom)
                if length is None or start < 0 or end > length or end <= start:
                    raise ValueError(f"Invalid region: {locus.locus}")
                values = np.asarray(signal.values(locus.chrom, start, end, numpy=True), dtype=float)
                if missing == "zero":
                    values = np.nan_to_num(values, nan=0.)
                edges = np.linspace(0, len(values), min(bins, len(values)) + 1, dtype=int)
                for left, right in zip(edges[:-1], edges[1:]):
                    segment = values[left:right]
                    finite = segment[np.isfinite(segment)]
                    records.append(dict(locus=locus.locus, track=track.track, scale_group=track.scale_group, chrom=locus.chrom,
                                        start=start + left, end=start + right, x=start + (left + right) / 2,
                                        value=finite.mean() if len(finite) else np.nan, covered_bases=len(finite)))
    return pd.DataFrame(records)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--tracks", type=Path, required=True)
    p.add_argument("--loci", type=Path, required=True)
    p.add_argument("--bins", type=int, default=450)
    p.add_argument("--missing", choices=["zero", "omit"], required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    tracks = read_table(a.tracks, ["track", "path", "scale_group"])
    tracks["path"] = [str((a.tracks.parent / value).resolve()) if not Path(value).is_absolute() else value for value in tracks.path]
    frame = locus_profiles(tracks, read_table(a.loci, ["locus", "chrom", "start", "end"]), a.bins, a.missing)
    out = new_output(a.output)
    save_table(frame, out / "tracks.tsv")
