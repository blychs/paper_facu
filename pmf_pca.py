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
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from load_data import load_data
from statsmodels.multivariate.pca import PCA

#%%
matrix, unc, meteo, gases, events, clusters = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                               'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                               'BA_events.xlsx', 'clusters.csv')

matrix = matrix[['Ag', 'Al', 'As', 'Ba', 'Ca', 'Cd', 'Cl', 'Co', 'Cr', 'Cu', 'EC',
       'Fe', 'K', 'Mg', 'Mn', 'Mo', 'NH4', 'NO3', 'Na no sol', 'Na sol', 'Ni', 'OC',
       'Pb', 'SO4', 'Sb', 'Se', 'Ti', 'V', 'Zn']]

#matrix = matrix[['Sb', 'Zn', 'Cu', 'Mo', 'Pb', 'Na sol', 'Cl']]
matrix = pd.concat([matrix, clusters], axis=1)


matrix = matrix.dropna()
#matrix.head()

pca_model = PCA(matrix.T)

fig, ax = plt.subplots(figsize=(15, 7))
lines = ax.plot(pca_model.factors.iloc[:, :3], lw=4, alpha=0.6)
ax.set_xticklabels(matrix.columns.values)
ax.axhline(0, color='k', linestyle='dashed')
#ax.set_xlabel("Element", size=17)
legend = fig.legend(lines, ["PC 1", "PC 2", "PC 3"], loc="center right")#, "PC 4", "PC 5", "PC 6"], loc="center right")
plt.show()

fig = pca_model.plot_scree(log_scale=False)
plt.show()

#%%

#%%
print(pca_model.factors.comp_09.sum())
#pca_model.loadings.plot.scatter(x="comp_00", y="comp_01")#, ax=ax)