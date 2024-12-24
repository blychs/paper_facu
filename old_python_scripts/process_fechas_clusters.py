# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: analysis
#     language: python
#     name: python3
# ---

# %%
import pandas as pd
import datetime as dt
import numpy as np

clusters = pd.read_csv("Fechas Clusters.txt", delim_whitespace=True)
#display(clusters)
clusters["YR"] = clusters["YR"] + 2000
clusters = clusters.rename(columns={"YR": "year",
                            "MO": "month", "DA":"day", "HR":"h"})
clusters["date"] = pd.to_datetime(clusters[["year", "month", "day", "h"]])
clusters["date"] = clusters["date"] - np.timedelta64(15, 'h')
clusters = clusters.sort_values("date", axis=0)
print(clusters)

#%%
length = len(clusters["date"].dt.date.unique())
dictval = {"date": sorted(clusters["date"].dt.date.unique()),
           "C1": np.zeros(length),
           "C2": np.zeros(length),
           "C3": np.zeros(length),
           "C4": np.zeros(length),
           "C5": np.zeros(length)}
outvalues = pd.DataFrame(dictval)

for date in clusters["date"].dt.date.unique():
    for j in clusters.loc[clusters["date"].dt.date == date]["CL#"]:
        if j == 1:
            outvalues.loc[outvalues["date"] == date, "C1"] += 1
        if j == 2:
            outvalues.loc[outvalues["date"] == date, "C2"] += 1
        if j == 3:
            outvalues.loc[outvalues["date"] == date, "C3"] += 1
        if j == 4:
            outvalues.loc[outvalues["date"] == date, "C4"] += 1
        if j == 5:
            outvalues.loc[outvalues["date"] == date, "C5"] += 1

print(outvalues)
outvalues.set_index(outvalues["date"], inplace=True)
outvalues.drop("date", axis=1, inplace=True)
outvalues.to_csv("clusters.csv")
print(clusters.loc[
    clusters['date'].dt.date == np.datetime64('2019-04-06')])