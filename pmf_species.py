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
#%%
print((matrix["NH4"]/matrix["PM2.5"] * 1e6).#/matrix["EC"]).
      #where(events["Event"].isin(["S", "SP", "SN"])).
      mean())

#print(events["Event"].isin(["S", "SN", "SP"]).sum())

#%%
for key in matrix[["SO4"]].keys():
    pollutant = key
    to_ug_g = lambda x: x/matrix["PM2.5"] * 1e6
    values = matrix.apply(to_ug_g)
    print(pollutant, f"mean concentration:\
           {(values[key]).mean():.1e}")
    print(confidence_interval(values.dropna(),
                              pollutant, 95,
                                n_resamples=int(1e5)).round(0))
    print("")

#%%
matrix["OC_EC"] = matrix["OC"]/matrix["EC"]
print(f"OC/EC\nmin = {matrix.OC_EC.min()}\n\
max = {matrix.OC_EC.max()}")
confidence_interval(matrix.
                    where(~events["Event"].isin(["S", "SP", "SN"])).
                    dropna(), "OC_EC", 95,
                     n_resamples=int(1e5))
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

meanprops = {"marker": "^",
             "markeredgecolor": "red",
             "markerfacecolor": "none",
             "markeredgewidth": 1,
            }
boxprops = {"facecolor": "none"}
medianprops = {"color": "tab:orange"}

PROPS = {"meanprops": meanprops,
         "boxprops": boxprops,
         "medianprops": medianprops}
        

with plt.style.context('seaborn-v0_8-paper'):
    over_gram = lambda x: x / (matrix["PM2.5"] / 1000000)
    fig,ax = plt.subplots(figsize=(8,5), ncols=5,
                          gridspec_kw={"width_ratios":[2,5,7,1,8]},
                          sharey=True)
    sns.boxplot(matrix[["OC", "EC"]].apply(over_gram),
        ax = ax[0], showmeans=True,
        **PROPS)
    ax[0].set_ylabel("µg/g")
    ax[0].set_yscale("log")
    ax[0].set_title("Carbonaceous")


    sns.boxplot(matrix[["SO4", "NO3", "Na$^+$", "Cl", "NH4"]].apply(over_gram),
                ax = ax[1], showmeans=True,
                **PROPS)
    ax[1].set_xticklabels(["SO$_4^{2-}$", "NO$_3^-$", "Na$^+$", "Cl$^-$",
                            "NH$_4^+$"])
    ax[1].set_title("SSA and II")

    sns.boxplot(matrix[["Ca", "Al", "Mg", "Fe", "Ti", "Mn"]].apply(over_gram),
                ax = ax[2], showmeans=True,
                **PROPS)
    ax[2].set_xticklabels(["Ca", "Al", "Mg", "Fe", "Ti", "Mn"])
    ax[2].set_title("Crustal (GM,\nconstruction/demolition)")

    sns.boxplot(matrix[["K"]].apply(over_gram),
                ax=ax[3], showmeans=True,
                **PROPS)
    ax[3].set_title("K")
    fig.tight_layout()
    plt.subplots_adjust(hspace=.0, wspace=.0)

    sns.boxplot(matrix[["Zn", "Ni", "Pb", "Cu", "Ba", "Sb", "V", "Mo"]].apply(over_gram),
        ax = ax[4], showmeans=True,
        **PROPS)
    ax[4].set_title("Traffic related and industrial\ntrace elements")

    fig.savefig("images/boxplot_by_source.png")
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