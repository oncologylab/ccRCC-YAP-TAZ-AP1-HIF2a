from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from figures._shared import DATA, read_table, numeric, bh, new_output, save_table


def downregulated(results, alpha, lfc, pcolumn):
    names = list(results)
    sets = {}
    universe = set()
    for name, frame in results.items():
        if frame.ENS_ID.duplicated().any():
            raise ValueError("Duplicate gene identifiers")
        universe.update(frame.ENS_ID)
        sets[name] = set(frame.loc[(frame[pcolumn] < alpha) & (frame.log2FoldChange < -lfc), "ENS_ID"])
    values = pd.DataFrame({"gene": sorted(universe)})
    for name, selected in sets.items():
        values[name] = values.gene.isin(selected)
    values["membership"] = ["+".join(name for name, flag in zip(names, flags) if flag) or "None"
                            for flags in values[names].itertuples(index=False, name=None)]
    counts = values.loc[values.membership != "None"].groupby("membership").size().rename("n_genes").reset_index()
    return values, counts


def jun_chromatin(frame):
    numeric(frame, ["control", "shJUNs"])
    if frame.duplicated(["stratum", "region"]).any():
        raise ValueError("Duplicate paired regions")
    records = []
    for stratum, part in frame.groupby("stratum"):
        test = wilcoxon(part.control, part.shJUNs, alternative="two-sided")
        records.append(dict(stratum=stratum, n_regions=len(part), ctrl_mean=part.control.mean(),
                            shJUNs_mean=part.shJUNs.mean(), ctrl_sd=part.control.std(), shJUNs_sd=part.shJUNs.std(), raw_p=test.pvalue))
    result = pd.DataFrame(records)
    result["FDR_p"] = bh(result.raw_p)
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("analysis", choices=["genes", "chromatin"])
    p.add_argument("--input", type=Path)
    p.add_argument("--data", type=Path, default=DATA)
    p.add_argument("--alpha", type=float)
    p.add_argument("--lfc", type=float)
    p.add_argument("--pcolumn", choices=["pvalue", "padj"])
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.analysis == "genes" and any(x is None for x in [a.alpha, a.lfc, a.pcolumn]):
        p.error("genes requires --alpha, --lfc and --pcolumn")
    if a.analysis == "chromatin" and a.input is None:
        p.error("chromatin requires --input")
    out = new_output(a.output)
    if a.analysis == "genes":
        results = {name: read_table(a.data / f"figure6d_DESeq2_{name}.tsv.gz", ["ENS_ID", "log2FoldChange", a.pcolumn])
                   for name in ["shYT", "shHIF2a", "shJUNs"]}
        values, counts = downregulated(results, a.alpha, a.lfc, a.pcolumn)
        save_table(values, out / "gene_membership.tsv")
        save_table(counts, out / "intersections.tsv")
    else:
        frame = read_table(a.input, ["stratum", "region", "control", "shJUNs"])
        save_table(jun_chromatin(frame), out / "tests.tsv")
        save_table(frame, out / "values.tsv")


if __name__ == "__main__":
    main()
