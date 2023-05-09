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
#import altair as alt
#from plotnine import ggplot, geom_smooth, aes, stat_smooth, facet_wrap
from load_data import load_data
import seaborn as sns
import matplotlib.pyplot as plt


matrix, unc, meteo, gases, events = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                               'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                               'BA_events.xlsx')


# %%
#alt.Chart(matrix.reset_index()).mark_line().encode(
#    x = 'date',
#    y = 'Ti'
#)

#%matplotlib widget

fig, ax = plt.subplots()#figsize=(7,7))
ax.errorbar(x = matrix.index, y = matrix["Ti"], yerr = unc["Ti"],
    capsize = 2, marker = '.', label = "Ti")
ax.errorbar(x = matrix.index, y = matrix["Sb"], yerr = unc["Sb"],
    capsize = 2, label = "Sb")
ax.legend()
fig.tight_layout()
plt.show()

over_gram = lambda x: x / (matrix["PM2.5"] / 1000000)
fig, ax = plt.subplots(figsize=(7.5,5), ncols=4)
sns.boxplot(matrix[["Ti", "Mg", "Al"]].apply(over_gram),
    ax = ax[0])
ax[0].set_ylabel("µg/g")
ax[0].set_yscale("log")
ax[0].set_title("Soil related")

sns.boxplot(matrix[["Sb", "Mo", "Cu"]].apply(over_gram),
    ax = ax[1])
ax[1].set_yscale("log")
ax[1].set_title("Break wear")

sns.boxplot(matrix[["Na sol", "Cl", "Mg"]].apply(over_gram),
            ax = ax[2])
ax[2].set_yscale("log")
ax[2].set_title("Sea Salt")

sns.boxplot(matrix[["OC", "EC"]].apply(over_gram),
            ax = ax[3])
ax[3].set_yscale("log")
ax[3].set_title("Carbonaceous")

fig.tight_layout()
plt.show()

#%%
fig = plt.figure()
ax = plt.axes()

ax.errorbar(
    x = matrix.index,
    y = matrix["Sb"],
    yerr = unc["Sb"],
    marker = '.',
    capsize = 2,
    capthick = 2,
)

print(matrix["Sb"].min())
print(matrix["Sb"].max())
print(matrix["Sb"].median())
print(matrix["Sb"].mean())
#%%
np.log10(matrix["Sb"]).plot.kde(label="Sb")
np.log10(matrix["Cu"]).plot.kde(label="Cu")
np.log10(matrix["PM2.5"]).plot.kde(label="PM$_{2.5}$")
np.log10(matrix["OC"]).plot.kde(label="OC")
np.log10(matrix["EC"]).plot.kde(label="EC")
np.log10(matrix["Na total"]).plot.kde(label="Na total")
plt.legend()