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
#     display_name: Python3 (analysis)
#     language: python
#     name: analysis
# ---

# %%
import pandas as pd
import numpy as np
import datetime as dt
from parserSFC import readSFC

data = pd.read_excel('PMF_BA_full.xlsx')
data['date'] = pd.to_datetime(data['date'])+ dt.timedelta(hours=12)
data.set_index(data['date'], inplace=True)
data.drop('date', axis=1, inplace=True)

meteo = readSFC('OBSERVATORIO.SFC')

orig_hours = data.index


for hour in orig_hours:
    for extended in range(1, 24):
        data.loc[hour + dt.timedelta(hours=extended)] = data.loc[hour]
data.sort_index(inplace=True)
meteo = meteo.loc[list(data.index)]

data = pd.concat([data, meteo[['ws', 'wd']]], axis=1)
data.to_csv('data_every_hour_duplicate.csv')


# %%
