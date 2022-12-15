# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# Load packages and data and convert negative values into **NaN**

# +
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
print(matrix.keys())

events = pd.read_excel('BA_events.xlsx', index_col='date')

# -

# Time series plot

# +
# %matplotlib widget
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
ax.set_xlabel("Date")
ax.set_ylabel("Mass concentration (µg/m$^3$)")
ax.legend()
plt.show()
# -



# +
# %matplotlib widget

methods = ['Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011']

#methods = ['Hand_2011']

d_methodQuality = {}


fig, axs = plt.subplots(3, 3, figsize=(20, 15))

i, j = 0, 0
for method in methods:
    d_methodQuality[method] = 0
    mass = mass_reconstruction(matrix, unc, equation=method)
    reconst = mass[0]/matrix['PM2.5']
    ureconst = np.sqrt((1/matrix['PM2.5'] * unc['PM2.5'])**2 + (matrix['PM2.5']/mass[0]/mass[0] * mass[2])**2)
    axs[i][j].errorbar(matrix.index, reconst,  yerr=ureconst, label=method, capsize=2, capthick=1)
    #axs[i][j].plot(matrix.index, reconst.where(np.logical_or(events['Event']=='S', events['Event']=='SP', events['Event']=='SN')), 'ro')
    axs[i][j].plot(matrix.index, reconst.where(events['Event']=='S'), 'ro', label='Smoke')
    axs[i][j].plot(matrix.index, reconst.where(events['Event']=='SP'), 'go', label='Smoke previous day')
    axs[i][j].plot(matrix.index, reconst.where(events['Event']=='SN'), 'ko', label='Smoke next day')
    axs[i][j].set_title(method)
    axs[i][j].legend()
    axs[i][j].axhline(0.8, color='k')
    axs[i][j].axhline(1.2, color='k')
    axs[i][j].tick_params(labelrotation=45)
    j += 1
    if j%3==0:
        j = 0
        i += 1
    d_methodQuality[method] = np.logical_and(((reconst + ureconst) > 0.8), ((reconst - ureconst) < 1.2)).sum()
    
plt.show()
    


# +
#Prepare correlation plots

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

        

