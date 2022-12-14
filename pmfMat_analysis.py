#!/usr/bin/env python
# coding: utf-8

# Load packages and data and convert negative values into **NaN**

# In[3]:


import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import pandas as pd
from funciones_pmfBA import mass_reconstruction

matrix = pd.read_excel('PMF_BA_full.xlsx', decimal=',', sheet_name='CONC', index_col='date')
matrix = matrix.rename(columns={'PM2,5': 'PM2.5'})
matrix[matrix < 0] = np.nan
matrix = matrix.reindex(sorted(matrix.columns), axis=1)

unc = pd.read_excel('PMF_BA_full.xlsx', decimal=',', sheet_name='UNC', index_col='date')
unc = unc.rename(columns={'PM2,5': 'PM2.5'})
unc[unc < 0] = np.nan
unc = unc.reindex(sorted(unc.columns), axis=1)


# Time series plot

# In[5]:


get_ipython().run_line_magic('matplotlib', 'widget')
mass = mass_reconstruction(matrix, unc, equation="Hand_2011")

plt.style.use('seaborn-v0_8-paper')

fig, ax = plt.subplots(figsize=(12,6))

#matrix['PM2.5'].plot(style='.-', label='PM2.5', ax=ax)
ax.errorbar(matrix.index, matrix['PM2.5'], yerr=unc['PM2.5'], marker='.', linestyle='-', capsize=3, capthick=1, label='PM2.5')
ax.errorbar(matrix.index, mass[1]['organic_mass'], yerr=mass[3]['uorganic_mass'], marker='.', capsize=3, capthick=1, linestyle='-', label="Organic mass")
ax.errorbar(matrix.index, mass[1]['inorganic_ions'], yerr=mass[3]['uinorganic_ions'], marker='.', capsize=3, capthick=1, linestyle='-', label="Inorganic ions")
ax.errorbar(matrix.index, mass[1]['geological_minerals'], yerr=mass[3]['ugeological_minerals'], marker='.', capsize=3, capthick=1, linestyle='-', label="Geological minerals")
ax.errorbar(matrix.index, mass[1]['elemental_C'], yerr=mass[3]['uelemental_C'], marker='.', capsize=3, capthick=1, linestyle='-', label="Elemental carbon")
ax.errorbar(matrix.index, mass[1]['salt'], yerr=mass[3]['usalt'], marker='.', capsize=3, capthick=1, linestyle='-', label="Sea salt")

ax.set_title('Time series')
ax.legend()
plt.show()


# In[5]:


(matrix.keys())


# In[8]:


for key1 in matrix.keys():
    i, j = 0, 0
    fig, axs = plt.subplots(6, 8, figsize=(45,30))
    for key2 in matrix.keys():
        axs[i][j].plot(matrix[key1], matrix[key2], 'o')
        axs[i][j].set_xlabel(key1)
        axs[i][j].set_ylabel(key2)
        j += 1
        if j%8 == 0:
            j=0
            i += 1
    fig.savefig(f'correlation_plots_{key1}.png')
    plt.close()

        

