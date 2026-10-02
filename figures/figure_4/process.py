from pathlib import Path
import argparse
import sys
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import read_table, numeric, new_output, save_table, read_bed, overlap_flags


def prism_values(path, table_id):
    root = ET.parse(path).getroot()
    for element in root.iter():
        element.tag = element.tag.rsplit("}", 1)[-1]
    table = root.find(f".//Table[@ID='{table_id}']")
    if table is None:
        raise ValueError(f"Missing Prism table {table_id}")
    columns = []
    for col in table.findall("YColumn"):
        values = [float(d.text) for d in col.findall("./Subcolumn/d")
                  if d.text and d.text.strip() and d.attrib.get("Excluded") != "1"]
        if values:
            columns.append(values)
    if len(columns) != 2:
        raise ValueError("Expected two populated Prism columns")
    return columns


def imaging(path):
    values, tests = [], []
    for table, panel, treated, method in [("Table0", "JUNB/HIF2a", "shYT", "exact-MW"), ("Table4", "JUNB/YAP", "shHIF2a", "Welch")]:
        a, b = prism_values(path, table)
        test = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact") if method == "exact-MW" else stats.ttest_ind(a, b, equal_var=False)
        tests.append(dict(panel=panel, n_ctrl=len(a), n_treated=len(b), method=method, pvalue=test.pvalue))
        for group, numbers in [("Ctrl", a), (treated, b)]:
            values.extend(dict(panel=panel, group=group, nucleus=i + 1, value=v) for i, v in enumerate(numbers))
    return pd.DataFrame(values), pd.DataFrame(tests)


def footprints(frame, flank):
    numeric(frame, ["position", "signal"])
    keys = ["panel", "condition", "site", "position"]
    if frame.duplicated(keys).any():
        raise ValueError("Duplicate site/position")
    frame = frame.loc[frame.position.abs() <= flank]
    if frame.empty:
        raise ValueError("No sites within flank")
    result = frame.groupby(["panel", "condition", "position"]).signal.agg(n_sites="count", mean="mean", sd="std", sem="sem").reset_index()
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("analysis", choices=["imaging", "footprints", "occupancy"])
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--yap", type=Path)
    p.add_argument("--hif2a", type=Path)
    p.add_argument("--jun", type=Path)
    p.add_argument("--flank", type=int, default=60)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.analysis == "occupancy" and not all([a.yap, a.hif2a, a.jun]):
        p.error("occupancy requires --yap, --hif2a and --jun")
    out = new_output(a.output)
    if a.analysis == "imaging":
        values, tests = imaging(a.input)
        save_table(values, out / "values.tsv")
        save_table(tests, out / "tests.tsv")
    elif a.analysis == "footprints":
        save_table(footprints(read_table(a.input, ["panel", "condition", "site", "position", "signal"]), a.flank), out / "footprints.tsv")
    else:
        regions = read_bed(a.input)
        for label, path in [("YAP", a.yap), ("HIF2a", a.hif2a), ("JUN", a.jun)]:
            regions[label] = overlap_flags(regions, read_bed(path))
        regions["class"] = ["+".join(name for name, flag in zip(["YAP", "HIF2a", "JUN"], flags) if flag) or "None"
                            for flags in regions[["YAP", "HIF2a", "JUN"]].itertuples(index=False, name=None)]
        save_table(regions, out / "occupancy.tsv")
        save_table(regions.groupby("class").size().rename("n_regions").reset_index(), out / "counts.tsv")


if __name__ == "__main__":
    main()
