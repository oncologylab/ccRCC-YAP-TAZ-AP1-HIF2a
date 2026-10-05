import numpy as np
import pandas as pd
from scipy.stats import chi2


def kaplan_meier(frame):
    if frame.duplicated(["group", "sample"]).any():
        raise ValueError("Duplicate sample")
    if not frame.event.isin([0, 1]).all() or (frame.time < 0).any():
        raise ValueError("Expected time >= 0 and event coded 0/1")
    if not np.isfinite(frame[["time", "event"]].to_numpy(float)).all():
        raise ValueError("Missing survival data")
    rows = []
    for group, part in frame.groupby("group", sort=False):
        survival = 1.0
        rows.append(dict(group=group, time=0., survival=survival, at_risk=len(part), events=0, censored=0))
        for time in np.sort(part.time.unique()):
            risk = int((part.time >= time).sum())
            events = int(((part.time == time) & (part.event == 1)).sum())
            censored = int(((part.time == time) & (part.event == 0)).sum())
            survival *= 1 - events / risk
            rows.append(dict(group=group, time=time, survival=survival, at_risk=risk, events=events, censored=censored))
    return pd.DataFrame(rows)


def logrank(frame):
    groups = frame.group.unique()
    if len(groups) != 2:
        raise ValueError("Two groups required")
    observed = expected = variance = 0.
    first = frame.group == groups[0]
    for time in np.sort(frame.loc[frame.event == 1, "time"].unique()):
        at_risk = frame.time >= time
        dead = (frame.time == time) & (frame.event == 1)
        n, n1, d, d1 = at_risk.sum(), (at_risk & first).sum(), dead.sum(), (dead & first).sum()
        observed += d1
        expected += d * n1 / n
        if n > 1:
            variance += n1 * (n - n1) * d * (n - d) / (n * n * (n - 1))
    statistic = (observed - expected) ** 2 / variance if variance > 0 else np.nan
    return pd.DataFrame([dict(group_a=groups[0], group_b=groups[1], chi2=statistic, pvalue=chi2.sf(statistic, 1))])
