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
from scipy.stats import lognorm
#import altair as alt
import plotnine as pn
from load_data import load_data
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import t #I'll use it to calculate the confidence interval


def describe_element(matrix, species):
    to_ug_g = lambda x: x/matrix["PM2.5"] * 1e6
    print(species, 'µg/g')
    print(matrix.apply(to_ug_g)[species].describe())


def confidence_interval(matrix, species, IC=95, n_resamples=10000):
    """Calculate the unbiased
    confidence interval
    based on
    https://towardsdatascience.com/confidence-intervals-with-python-bfa28ebb81c
    """
    values = matrix[species].values
    confidence = IC/100
    bootstrap_values = [np.random.choice(values, size=int(len(values)), replace=True).mean()
                        for i in range(n_resamples)]
    ic = np.percentile(bootstrap_values, [100*(1-confidence)/2, 100*(1-(1-confidence)/2)])
    return(ic)
    

matrix, unc, meteo, gases, events, clusters = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                               'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                               'BA_events.xlsx', 'clusters.csv')


matrix.rename(columns={'Na sol': 'Na$^+$'}, inplace=True)
#clusters = clusters.astype(int)
unc = unc.add_prefix('unc_')
matrix_withunc = pd.concat([matrix, unc, events, clusters], axis=1)
matrix_withunc["C1"] = matrix_withunc["C1"].astype(float)
#print(matrix_withunc)

pollutant = "OC"
#print(matrix["OC"].describe())
to_ug_g = lambda x: x/matrix["PM2.5"] * 1e6
#print(matrix.apply(to_ug_g)[pollutant])
print(confidence_interval(matrix.apply(to_ug_g).dropna(), pollutant, 95, n_resamples=int(1e5)))
#display(matrix_withunc)


#%%
print(list(matrix.keys()))
with plt.style.context('ggplot'):
    fig, ax = plt.subplots()#figsize=(7,7))
    ax.errorbar(x = matrix.index, y = 15* matrix["Ba"], yerr = unc["unc_Ti"],
        capsize = 2, marker = '.', label = "Ba")
    ax.errorbar(x = matrix.index, y = 0.5*matrix["Ca"], yerr = unc["unc_Sb"],
        capsize = 2, label = "Ca")
    ax.errorbar(x = matrix.index, y = 0.05 * matrix["Na total"], yerr = unc["unc_Sb"],
        capsize = 2, label = "Na tot")
    ax.errorbar(x = matrix.index, y = matrix["Mg"], yerr = unc["unc_Sb"],
        capsize = 2, label = "Mg")
    ax.legend()
    fig.tight_layout()
    plt.show()

#%%
with plt.style.context('seaborn-v0_8-paper'):
    over_gram = lambda x: x / (matrix["PM2.5"] / 1000000)
    fig,ax = plt.subplots(figsize=(7.5,5), ncols=5,
                          gridspec_kw={"width_ratios":[3,6,2,2,2]},
                          sharey=True)
    sns.boxplot(matrix[["Ti", "Mg", "Al"]].apply(over_gram),
        ax = ax[0])
    ax[0].set_ylabel("µg/g")
    ax[0].set_yscale("log")
    ax[0].set_title("Crustal")

    sns.boxplot(matrix[["Sb", "Mo", "Cu", "Pb", "V", "Zn"]].apply(over_gram),
        ax = ax[1])
    ax[1].set_title("Traffic related")

    sns.boxplot(matrix[["Na$^+$", "Cl"]].apply(over_gram),
                ax = ax[2])
    ax[2].set_xticklabels(["Na$^+$", "Cl$^-$"])
    ax[2].set_title("Sea Salt")

    sns.boxplot(matrix[["OC", "EC"]].apply(over_gram),
                ax = ax[3])
    ax[3].set_title("Carbonaceous")

    sns.boxplot(matrix[["K", "Na no sol"]].apply(over_gram),
                ax=ax[4])
    ax[4].set_title("Other\nmarkers")
    fig.tight_layout()
    plt.subplots_adjust(hspace=.0, wspace=.0)
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
plt.show()
#%%
np.log10(matrix["Sb"]).plot.kde(label="Sb")
np.log10(matrix["Cu"]).plot.kde(label="Cu")
np.log10(matrix["PM2.5"]).plot.kde(label="PM$_{2.5}$")
np.log10(matrix["OC"]).plot.kde(label="OC")
np.log10(matrix["EC"]).plot.kde(label="EC")
np.log10(matrix["Na total"]).plot.kde(label="Na total")
plt.xticks()
plt.legend()
plt.show()

#%%
import numpy as np, seaborn as sns, matplotlib.pyplot as plt
np.random.seed(1)
data = np.power(np.random.randn(1000), 10)

fig, ax = plt.subplots()
sns.kdeplot(np.log10(matrix["Sb"]), fill=True, ax=ax)
fig.canvas.draw()
locs, labels = plt.xticks()
# u2212 is the matplotlib's medium dash for negative numbers.
ax.set(xticklabels=[10 ** int(i.get_text().replace(u'\u2212', '-'))
                    for i in labels])
# Or for scientific notation:
# ax.set(xticklabels=["$10^" + i.get_text() + "$" for i in labels])
plt.show()

#%%
fig, ax = plt.subplots()
ax.errorbar(x=matrix.index, y=matrix["Sb"], yerr=unc["unc_Sb"],
            marker='.', capsize=2, label="Sb")
ax.errorbar(x=matrix.index, y=matrix["Cu"], yerr=unc["unc_Cu"],
            marker='.', capsize=2, label="Cu")
ax.legend()

plt.show()
#%%
#%matplotlib widget
print(matrix_withunc.loc[
    (matrix.index.date > np.datetime64('2019-05-10')) & 
    (matrix.index.date <= np.datetime64('2019-07-01'))][['VentCoef', 'ws', 'C1', 'C2', 'C3', 'C4', 'C5']])
#%%
fig, ax = plt.subplots()

matrix["Cl"].plot(label='Ca')
matrix["Na total"].plot(label="Na tot")
#(0.0004 * matrix['VentCoef']).plot(label='VentCoef')
#(0.1 * matrix["PM2.5"]).plot(label="PM2.5")
plt.legend()