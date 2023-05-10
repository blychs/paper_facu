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
#import plotnine as pn
from load_data import load_data
import seaborn as sns
import matplotlib.pyplot as plt


matrix, unc, meteo, gases, events = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                               'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                               'BA_events.xlsx')


unc = unc.add_prefix('unc_')
matrix_withunc = pd.concat([matrix, unc, events], axis=1)
#display(matrix_withunc)
# %%
#alt.Chart(matrix.reset_index()).mark_line().encode(
#    x = 'date',
#    y = 'Ti'
#)

#%matplotlib widget

p = (
    pn.ggplot(matrix_withunc.reset_index().rename(columns={"PM2.5": "PM25",
                                                           "unc_PM2.5": "unc_PM25"})) +
    #pn.geom_line(pn.aes(x = "date", y = "Ti"), color='b') +
    #pn.geom_errorbar(pn.aes(x = "date", ymin = "Ti-unc_Ti",
    #                         ymax = "Ti+unc_Ti", width=4), color='k') +
    #pn.geom_point(pn.aes(x = "date", y = "Ti"), color='b') +
    pn.geom_line(pn.aes(x = "date", y = "PM25"), color='k') +
    pn.geom_point(pn.aes(x = "date", y = "PM25", color="Event")) +
    pn.geom_errorbar(pn.aes(x = "date", ymin = f"PM25-unc_PM25",
                             ymax = f"PM25+unc_PM25", width=4, color="Event")) +
    pn.scale_y_continuous(trans="log10")
)

p.draw()
plt.show()

with plt.style.context('ggplot'):
    fig, ax = plt.subplots()
    ax.plot(matrix.index, matrix["PM2.5"], '.-r')


#%%
print(list(matrix.keys()))
with plt.style.context('ggplot'):
    fig, ax = plt.subplots()#figsize=(7,7))
    ax.errorbar(x = matrix.index, y = matrix["Ti"], yerr = unc["unc_Ti"],
        capsize = 2, marker = '.', label = "Ti")
    ax.errorbar(x = matrix.index, y = matrix["Sb"], yerr = unc["unc_Sb"],
        capsize = 2, label = "Sb")
    ax.legend()
    fig.tight_layout()
    plt.show()

#%%
with plt.style.context('seaborn-v0_8-paper'):
    over_gram = lambda x: x / (matrix["PM2.5"] / 1000000)
    fig,ax = plt.subplots(figsize=(7.5,5), ncols=5,
                           gridspec_kw={"width_ratios":[3,5,3,2,2]})
    sns.boxplot(matrix[["Ti", "Mg", "Al"]].apply(over_gram),
        ax = ax[0])
    ax[0].set_ylabel("µg/g")
    ax[0].set_yscale("log")
    ax[0].set_title("Soil related")

    sns.boxplot(matrix[["Sb", "Mo", "Cu", "Pb", "Zn"]].apply(over_gram),
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

    sns.boxplot(matrix[["V", "Zn"]], ax=ax[4])
    ax[4].set_yscale("log")
    ax[4].set_title("Exhaust")
    fig.tight_layout()
    plt.show()

#%%
fig, ax = plt.subplots()

ax.errorbar(
    x = matrix.index,
    y = matrix["Sb"],
    yerr = unc["unc_Sb"],
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