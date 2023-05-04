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
import numpy as np
import datetime as dt
import pandas as pd
import altair as alt
from load_data import load_data
import seaborn as sns


matrix, unc, meteo, gases, events = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                               'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                               'BA_events.xlsx')

# %%
alt.Chart(matrix.reset_index()).mark_line().encode(
    x = 'date',
    y = 'Ti'
)

#%%
line = alt.Chart(matrix.reset_index()).mark_line().encode(
    x = 'date',
    y = 'Ti'
)

ebars = alt.Chart(unc.reset_index()).mark_errorbar().encode(
    x = 'date',
    y = 'Ti'
)

line + ebars

#%%
