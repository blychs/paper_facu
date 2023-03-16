# -*- coding: utf-8 -*-
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

# %% [markdown]
# Load packages and data and convert negative values into **NaN**

# %%
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
from scipy.stats import linregress, spearmanr, zscore
import statsmodels.api as sm
from funciones_pmfBA import mass_reconstruction, mass_reconstruction_mod, percentage_with_err
from funciones_pmfBA import estimation_om_oc, calculate_seasonal
from load_data import load_data


matrix, unc, meteo, gases, events = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                               'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                               'BA_events.xlsx')
display(matrix)


# %%
pd.set_option('display.float_format', '{:.2g}'.format)

matrix_seasonal = calculate_seasonal(matrix)

#print(matrix_seasonal.to_latex())


# %%
# %matplotlib widget
fig, ax = plt.subplots()

ax.errorbar(matrix.index, matrix['PM2.5'], yerr=unc['PM2.5'], capsize=2, capthick=1, marker='.', zorder=0)
ax.plot(matrix.index, matrix['PM2.5'].where(events['Event'].isin(['S'])),  'o', label='Smoke')
ax.plot(matrix.index, matrix['PM2.5'].where(events['Event'].isin(['SP'])), 'o', label='Smoke previous day')
ax.plot(matrix.index, matrix['PM2.5'].where(events['Event'].isin(['SN'])), 'o', label='Smoke next day')
ax.axhline(50,   color='c', linestyle='-.', label='Interim target 2')
ax.axhline(37.5, color='g', linestyle='--', label='Interim target 3')
ax.axhline(35,   color='y', linestyle=':',  label='Local regulations')
ax.axhline(25,   color='m', linestyle='-.', label='Interim target 4')
ax.axhline(15,   color='r', linestyle='--', label='AQ Guideline')
ax.legend()
ax.set_xlabel('Date')
ax.set_ylabel('PM$_{2.5}$ mass concentration (µg m$^{-3}$)')
# ax.grid()
fig.savefig('HV_PM25.png')
plt.show()

## Calculate num exceedances (percentage)
valid_filters = matrix['PM2.5'].notna().sum()
print(valid_filters)
exceedance = lambda target: ((matrix['PM2.5'] > target).sum()/valid_filters * 100).round(2)

print('Interim target 1:', exceedance(75))
print('Interim target 2:', exceedance(50))
print('Interim target 3:', exceedance(37.5))
print('Interim target 4:', exceedance(25))
print('Interim target AQG:', exceedance(15))
print('Local regulation (35 µg/m$^3$) =', exceedance(35))
print('Yearly mean:', matrix['PM2.5'].mean().round(2))

# %%
matrix['month'] = pd.DatetimeIndex(matrix.index)
matrix['month'] = matrix['month'].dt.to_period('M')


#print(matrix.groupby(matrix['month']).mean().round(2)[['PM2.5', 'ws', 'VentCoef', 'temp']].style.to_latex())

fig, ax = plt.subplots()

### Boxplots, width of box depends on number of data points.

bins, groups = zip(*matrix['PM2.5'].dropna().groupby(matrix['month']))
lengths = np.array([len(group) for group in groups])
print(bins)
max_width = 0.8
ax.boxplot(groups, widths=max_width * lengths / lengths.max(),
            patch_artist=True, boxprops={'fill': None})
ax.set_ylabel('PM$_{2.5}$ concentration (µg m$^{-3}$)')
ax.set_xlabel('Month')
ax.set_xticklabels(bins, rotation=45, ha='right')
ax.grid()
fig.tight_layout()
fig.savefig('PM_boxplot.png')
plt.show()

# %%
#print(plt.rcParams.keys())

fig, ax = plt.subplots()

ax.scatter(matrix['VentCoef'], matrix['PM2.5'], label='No event')
ax.scatter(matrix['VentCoef'], matrix['PM2.5'].where(events['Event'].isin(['S', 'SP', 'SN'])), label='Event')
ax.set_xlabel('Ventilation Coefficient (m$^2$/s)')
ax.set_ylabel('PM2.5 (µg/m$^3$)')
ax.legend()
fig.savefig('PM2.5_vs_vent.png')
plt.show()

#spear


with plt.rc_context({'axes.labelsize': 15}):
    spearman_corr = matrix.corr(numeric_only=True, method='spearman')
    plt.figure(figsize=(35,20))
    plt.tick_params(axis='x', which='both', labelsize=10, labelrotation=90, labelbottom = True, bottom=True, top = True, labeltop=True)
    plt.tick_params(axis='y', which='both', labelsize=10, labelleft = True, left=True, labelright = True, right=True)
    sns.heatmap(spearman_corr, annot=True,  cmap='RdBu_r', vmin=-1, vmax=1)
    plt.title('Spearman correlation', fontsize=20)
    plt.savefig('heatmap_spearman_corr.png')
    plt.show()

# %%

with plt.rc_context({'axes.labelsize': 15}):
    spearman_corr = matrix.corr(numeric_only=True, method='pearson')
    plt.figure(figsize=(35,20))
    plt.title('Pearson correlation', fontsize=20)
    plt.tick_params(axis='x', which='both', labelsize=10, labelrotation=90, labelbottom = True, bottom=True, top = True, labeltop=True)
    plt.tick_params(axis='y', which='both', labelsize=10, labelleft = True, left=True, labelright = True, right=True)
    sns.heatmap(spearman_corr, annot=True,  cmap='RdBu_r', vmin=-1, vmax=1)
    plt.savefig('heatmap_pearson_corr.png')
    plt.show()

# %%
fig, ax = plt.subplots()
ax.scatter(events['AOD440'], events['Alpha'])
ax.scatter(events['AOD440'], events['Alpha'].where(events['Event'].isin(['S', 'SP', 'SN'])), label='Event')
ax.axhline(0.85, color='k', linestyle='dashed')
ax.axvline(0.185, color='k', linestyle='dashed')
ax.set_xlabel('AOD 440nm', size=10)
ax.set_ylabel(r'$\alpha$', size=10)
ax.legend(loc=4)
fig.savefig('AOD_alpha.png')
plt.show()

# %%
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

# %%
# %matplotlib widget
mass = mass_reconstruction(matrix, unc, equation="Hand_2011")

organic_mass_per = percentage_with_err(mass[1]['organic_mass'], mass[0], mass[3]['uorganic_mass'], mass[2])
inorganic_ions_per = percentage_with_err(mass[1]['inorganic_ions'], mass[0], mass[3]['uinorganic_ions'], mass[2])
geological_minerals_per = percentage_with_err(mass[1]['geological_minerals'], mass[0], mass[3]['ugeological_minerals'], mass[2])
EC_per = percentage_with_err(mass[1]['elemental_C'], mass[0], mass[3]['uelemental_C'], mass[2])
ssa_per = percentage_with_err(mass[1]['salt'], mass[0], mass[3]['usalt'], mass[2])

plt.style.use('seaborn-v0_8-paper')

#fig, ax = plt.subplots(figsize=(12,6))
#
#
##ax.errorbar(matrix.index, matrix['PM2.5'], yerr=unc['PM2.5'], marker='.', linestyle='-', capsize=3, capthick=1, label='PM2.5')
#ax.errorbar(matrix.index, organic_mass_per['perc'], yerr=organic_mass_per['uperc'], marker='.', capsize=3, capthick=1, linestyle='-', label="Organic mass")
#ax.errorbar(matrix.index, inorganic_ions_per['perc'], yerr=inorganic_ions_per['uperc'], marker='.', capsize=3, capthick=1, linestyle='-', label="Inorganic ions")
#ax.errorbar(matrix.index, geological_minerals_per['perc'], yerr=geological_minerals_per['uperc'], marker='.', capsize=3, capthick=1, linestyle='-', label="Geological minerals")
#ax.errorbar(matrix.index, EC_per['perc'], yerr=EC_per['uperc'], marker='.', capsize=3, capthick=1, linestyle='-', label="Elemental carbon")
#ax.errorbar(matrix.index, mass[1]['salt'], yerr=mass[3]['usalt'], marker='.', capsize=3, capthick=1, linestyle='-', label="Sea salt")
#
#ax.set_title('Time series as percentage of reconstructed mass')
#ax.set_xlabel("Date")
#ax.set_ylabel("Mass concentration (µg/m$^3$)")
#ax.legend()
#plt.show()

reconst = mass[0]/matrix['PM2.5'] * 100
ureconst = np.sqrt((1/matrix['PM2.5'] * unc['PM2.5'])**2 + (matrix['PM2.5']/mass[0]/mass[0] * mass[2])**2) * 100
smoke_dates = list(matrix.index.where(events['Event']=='S').dropna())
print(smoke_dates)


width=2.5

fig, ax = plt.subplots(nrows=4, ncols=2, figsize=(15, 12), sharex=True)

plt.subplots_adjust(hspace=.0)
plt.subplots_adjust(wspace=.0)

fig.suptitle('Mass reconstruction')
ax[0][0].set_title('Hand 2011')
ax[0][0].bar(matrix.index, organic_mass_per['perc'], width, yerr=organic_mass_per['uperc'], label='Organic mass')
ax[0][0].bar(matrix.index, inorganic_ions_per['perc'], width, yerr=inorganic_ions_per['uperc'], bottom=organic_mass_per['perc'],label='Inorganic ions')
ax[0][0].bar(matrix.index, geological_minerals_per['perc'], width, yerr=geological_minerals_per['uperc'],
         bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc']),label='Geological minerals')
ax[0][0].bar(matrix.index, EC_per['perc'], width, yerr=EC_per['uperc'],
          bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + geological_minerals_per['perc']),label='EC')
ax[0][0].bar(matrix.index, ssa_per['perc'], width, yerr=ssa_per['uperc'],
          bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + geological_minerals_per['perc'] + EC_per['perc']),label='SSA')
ax[0][0].set_ylabel('Category of rec. mass (%)')
ax[0][0].legend(loc=4)

ax[1][0].bar(matrix.index, mass[1]['organic_mass'], width, yerr=mass[3]['uorganic_mass'], label='Organic mass')
ax[1][0].bar(matrix.index, mass[1]['inorganic_ions'], width, yerr=mass[3]['uinorganic_ions'],
             bottom=mass[1]['organic_mass'],label='Inorganic ions')
ax[1][0].bar(matrix.index, mass[1]['geological_minerals'], width, yerr=mass[3]['ugeological_minerals'],
         bottom=(mass[1]['organic_mass'] + mass[1]['inorganic_ions']),label='Geological minerals')
ax[1][0].bar(matrix.index, mass[1]['elemental_C'], width, yerr=mass[3]['uelemental_C'],
          bottom=(mass[1]['organic_mass'] + mass[1]['inorganic_ions'] + mass[1]['geological_minerals']),label='EC')
ax[1][0].bar(matrix.index, mass[1]['salt'], width, yerr=mass[3]['usalt'],
          bottom=(mass[1]['organic_mass'] + mass[1]['inorganic_ions'] + mass[1]['geological_minerals'] + mass[1]['elemental_C']),label='SSA')
ax[1][0].errorbar(matrix.index, matrix['PM2.5'], yerr=unc['PM2.5'], capsize=2, capthick=1, marker='.', ecolor='tan', color='tan')
ax[1][0].set_ylabel('Category of rec. mass ($\mu$g/m$^3$)')
ax[1][0].legend(loc=1)


ax[2][0].errorbar(x=matrix.index, y=reconst,  yerr=ureconst, capsize=2, capthick=1, marker='.', ecolor='cornflowerblue')
ax[2][0].axhspan(80, 120, facecolor='y', alpha=0.5)
ax[2][0].plot(matrix.index, reconst.where(events['Event']=='S'), 'o', label='Smoke')
ax[2][0].plot(matrix.index, reconst.where(events['Event']=='SP'), 'o', label='Smoke prev. day')
ax[2][0].plot(matrix.index, reconst.where(events['Event']=='SN'), 'o', label='Smoke next day')
ax[2][0].axhline(100, color='k')
ax[2][0].legend()
ax[2][0].set_ylabel('Percentage of PM2.5 reconstructed')

ax[3][0].errorbar(matrix.index, matrix['PM2.5'], yerr=unc['PM2.5'], capsize=2, capthick=1, marker='.', label='Measured PM2.5')
ax[3][0].errorbar(matrix.index, mass[0], mass[2], capsize=2, capthick=1, marker='.', label='Reconstructed PM2.5') 
ax[3][0].legend()
ax[3][0].set_ylabel('Mass $\mu$g/m$^3$')


mass = mass_reconstruction_mod(matrix, unc, events=events, equation="DeBell_2006")
reconst = mass[0]/matrix['PM2.5'] * 100
ureconst = np.sqrt((1/matrix['PM2.5'] * unc['PM2.5'])**2 + (matrix['PM2.5']/mass[0]/mass[0] * mass[2])**2) * 100
organic_mass_per = percentage_with_err(mass[1]['organic_mass'], mass[0], mass[3]['uorganic_mass'], mass[2])
inorganic_ions_per = percentage_with_err(mass[1]['inorganic_ions'], mass[0], mass[3]['uinorganic_ions'], mass[2])
geological_minerals_per = percentage_with_err(mass[1]['geological_minerals'], mass[0], mass[3]['ugeological_minerals'], mass[2])
EC_per = percentage_with_err(mass[1]['elemental_C'], mass[0], mass[3]['uelemental_C'], mass[2])
#ssa_per = percentage_with_err(mass[1]['salt'], mass[0], mass[3]['usalt'], mass[2])

ax[0][1].set_title('DeBell 2006')
ax[0][1].bar(matrix.index, organic_mass_per['perc'], width, yerr=organic_mass_per['uperc'], label='Organic mass')
ax[0][1].bar(matrix.index, inorganic_ions_per['perc'], width, yerr=inorganic_ions_per['uperc'], bottom=organic_mass_per['perc'],label='Inorganic ions')
ax[0][1].bar(matrix.index, geological_minerals_per['perc'], width, yerr=geological_minerals_per['uperc'],
          bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc']),label='Geological minerals')
ax[0][1].bar(matrix.index, EC_per['perc'], width, yerr=EC_per['uperc'],
          bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + geological_minerals_per['perc']),label='EC')
#ax[0][1].bar(matrix.index, ssa_per['perc'], width, yerr=ssa_per['uperc'],
#          bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + geological_minerals_per['perc'] + EC_per['perc']),label='SSA')
ax[0][1].yaxis.tick_right() 
ax[0][1].yaxis.set_label_position("right")
ax[0][1].set_ylabel('Category of rec. mass (%)')
ax[0][1].legend(loc=4)

ax[1][1].bar(matrix.index, mass[1]['organic_mass'], width, yerr=mass[3]['uorganic_mass'], label='Organic mass')
ax[1][1].bar(matrix.index, mass[1]['inorganic_ions'], width, yerr=mass[3]['uinorganic_ions'],
             bottom=mass[1]['organic_mass'],label='Inorganic ions')
ax[1][1].bar(matrix.index, mass[1]['geological_minerals'], width, yerr=mass[3]['ugeological_minerals'],
         bottom=(mass[1]['organic_mass'] + mass[1]['inorganic_ions']),label='Geological minerals')
ax[1][1].bar(matrix.index, mass[1]['elemental_C'], width, yerr=mass[3]['uelemental_C'],
          bottom=(mass[1]['organic_mass'] + mass[1]['inorganic_ions'] + mass[1]['geological_minerals']),label='EC')
#ax[1][1].bar(matrix.index, mass[1]['salt'], width, yerr=mass[3]['usalt'],
#          bottom=(mass[1]['organic_mass'] + mass[1]['inorganic_ions'] + mass[1]['geological_minerals'] + mass[1]['elemental_C']),label='SSA')
ax[1][1].errorbar(matrix.index, matrix['PM2.5'], yerr=unc['PM2.5'], capsize=2, capthick=1, marker='.', ecolor='tan', color='tan')
ax[1][1].yaxis.tick_right() 
ax[1][1].yaxis.set_label_position("right")
ax[1][1].set_ylabel('Category of rec. mass ($\mu$g/m$^3$)')
ax[1][1].legend(loc=1)

ax[2][1].errorbar(x=matrix.index, y=reconst,  yerr=ureconst, capsize=2, capthick=1, marker='.', ecolor='cornflowerblue')
ax[2][1].axhspan(80, 120, facecolor='y', alpha=0.5)
ax[2][1].plot(matrix.index, reconst.where(events['Event']=='S'), 'o', label='Smoke')
ax[2][1].plot(matrix.index, reconst.where(events['Event']=='SP'), 'o', label='Smoke prev. day')
ax[2][1].plot(matrix.index, reconst.where(events['Event']=='SN'), 'o', label='Smoke next day')
ax[2][1].axhline(100, color='k')
ax[2][1].legend(loc=1)
ax[2][1].yaxis.tick_right() 
ax[2][1].yaxis.set_label_position("right")
ax[2][1].set_ylabel('Percentage of PM2.5 reconstructed')

ax[3][1].errorbar(matrix.index, matrix['PM2.5'], yerr=unc['PM2.5'], capsize=2, capthick=1, marker='.', label='Measured PM2.5')
ax[3][1].errorbar(matrix.index, mass[0], mass[2], capsize=2, capthick=1, marker='.', label='Reconstructed PM2.5') 
ax[3][1].legend()
ax[3][1].yaxis.tick_right() 
ax[3][1].yaxis.set_label_position("right")
ax[3][1].set_ylabel('Mass $\mu$g/m$^3$')

fig.savefig('stacked_massreconst.png')


plt.show()

# %%
# %matplotlib widget
mass = mass_reconstruction(matrix, unc, equation="Hand_2011")

values = mass[1]['inorganic_ions'] + mass[1]['geological_minerals'] + mass[1]['elemental_C'] + mass[1]['salt']
uvalues = mass[3]['uinorganic_ions'] + mass[3]['ugeological_minerals'] + mass[3]['uelemental_C'] + mass[3]['usalt'] 


plt.style.use('seaborn-v0_8-paper')

fig, ax = plt.subplots(figsize=(15,7.5))

#matrix['PM2.5'].plot(style='.-', label='PM2.5', ax=ax)
ax.errorbar(matrix.index, matrix['PM2.5'], yerr=unc['PM2.5'], marker='.', linestyle='-', capsize=3, capthick=1, label='PM2.5')
ax.errorbar(matrix.index, mass[1]['organic_mass'], yerr=mass[3]['uorganic_mass'], marker='.', capsize=3, capthick=1, linestyle='-', label="Organic mass")
ax.errorbar(matrix.index, values, yerr=uvalues, marker='.', capsize=3, capthick=1, linestyle='-', label="Other")
ax.errorbar(matrix.index, mass[0], yerr=mass[2], marker='.', capsize=3, capthick=1, linestyle='-', label="Closure")

ax.set_title('Time series')
ax.set_xlabel("Date")
ax.set_ylabel("Mass concentration (µg/m$^3$)")
ax.legend()
plt.show()

# %%
# #%matplotlib widget

methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011', 'Simon_2011']

#methods = ['Hand_2011']

d_methodQuality = {}

plt.style.use('seaborn-v0_8-paper')
fig, axs = plt.subplots(4, 3, figsize=(16, 12))

i, j = 0, 0
for method in methods:
    d_methodQuality[method] = 0
    mass = mass_reconstruction(matrix, unc, equation=method)
    reconst = mass[0]/matrix['PM2.5'] * 100
    ureconst = np.sqrt((1/matrix['PM2.5'] * unc['PM2.5'])**2 + (matrix['PM2.5']/mass[0]/mass[0] * mass[2])**2) * 100
    axs[i][j].errorbar(matrix.index, reconst,  yerr=ureconst, capsize=2, capthick=1, marker='.', ecolor='cornflowerblue', zorder=0)
    axs[i][j].plot(matrix.index, reconst.where(events['Event']=='S'), 'o', label='Smoke', zorder=1)
    axs[i][j].plot(matrix.index, reconst.where(events['Event']=='SP'), 'D', label='Smoke previous day', zorder=2)
    axs[i][j].plot(matrix.index, reconst.where(events['Event']=='SN'), 's', label='Smoke next day', zorder=3)
    axs[i][j].set_title(method)
    axs[i][j].legend(loc=9)
    axs[i][j].axhline(80, color='k')
    axs[i][j].axhline(120, color='k')
    axs[i][j].tick_params(labelrotation=0)
    j += 1
    if j%3==0:
        j = 0
        i += 1
    d_methodQuality[method] = np.logical_and(((reconst + ureconst) > 80), ((reconst - ureconst) < 120)).sum()

for x in range(0, 3):
    axs[x][0].set_ylabel('Mass reconstructed [%]')
    axs[-1][x].set_xlabel('Date')
    
plt.show()
print(d_methodQuality)


# %%
# #%matplotlib widget

methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011', 'Simon_2011']

#methods = ['Hand_2011']

d_methodQuality = {}

plt.style.use('seaborn-v0_8-paper')
fig, axs = plt.subplots(4, 3, figsize=(16, 12))

i, j = 0, 0
for method in methods:
    d_methodQuality[method] = 0
    mass = mass_reconstruction(matrix, unc, equation=method)
    reconst = mass[0]/matrix['PM2.5'] * 100
    ureconst = np.sqrt((1/matrix['PM2.5'] * unc['PM2.5'])**2 + (matrix['PM2.5']/mass[0]/mass[0] * mass[2])**2) * 100
    axs[i][j].errorbar(matrix.index, reconst,  yerr=ureconst, capsize=2, capthick=1, marker='.', ecolor='cornflowerblue', zorder=0)
#    axs[i][j].plot(matrix.index, reconst.where(events['AOD440']>0.185), 'or', label='AOD', zorder=1)
#    axs[i][j].plot(matrix.index, reconst.where(events['Alpha']>0.85), 'Xg', label='Alpha', zorder=3)
    axs[i][j].plot(matrix.index, reconst.where(events['Event']=="S"), 's', label='Smoke', zorder=1)
    axs[i][j].plot(matrix.index, reconst.where(events['Alpha']>0.85).where(events['AOD440']>0.185), '^g', label='Alpha+AOD', zorder=2)
    axs[i][j].set_title(method)
    axs[i][j].legend(loc=9)
    axs[i][j].axhline(80, color='k')
    axs[i][j].axhline(120, color='k')
    axs[i][j].tick_params(labelrotation=0)
    j += 1
    if j%3==0:
        j = 0
        i += 1
    d_methodQuality[method] = np.logical_and(((reconst + ureconst) > 80), ((reconst - ureconst) < 120)).sum()

for x in range(0, 3):
    axs[x][0].set_ylabel('Mass reconstructed [%]')
    axs[-1][x].set_xlabel('Date')
    
plt.show()
print(d_methodQuality)


# %%
# #%matplotlib widget

methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011', 'Simon_2011']

#methods = ['Hand_2011']

d_methodQuality = {}

plt.style.use('seaborn-v0_8-paper')
fig, axs = plt.subplots(4, 3, figsize=(16, 12))

i, j = 0, 0
for method in methods:
    d_methodQuality[method] = 0
    mass = mass_reconstruction_mod(matrix, unc, events, equation=method)
    reconst = mass[0]/matrix['PM2.5'] * 100
    ureconst = np.sqrt((1/matrix['PM2.5'] * unc['PM2.5'])**2 + (matrix['PM2.5']/mass[0]/mass[0] * mass[2])**2) * 100
    axs[i][j].errorbar(matrix.index, reconst,  yerr=ureconst, capsize=2, capthick=1, marker='.', ecolor='cornflowerblue', zorder=0)
#    axs[i][j].plot(matrix.index, reconst.where(events['AOD440']>0.185), 'or', label='AOD', zorder=1)
#    axs[i][j].plot(matrix.index, reconst.where(events['Alpha']>0.85), 'Xg', label='Alpha', zorder=3)
    axs[i][j].plot(matrix.index, reconst.where(events['Event']=="S"), 's', label='Smoke', zorder=1)
    axs[i][j].plot(matrix.index, reconst.where(events['Alpha']>0.85).where(events['AOD440']>0.185), '^g', label='Alpha+AOD', zorder=2)
    axs[i][j].set_title(method)
    axs[i][j].legend(loc=9)
    axs[i][j].axhline(80, color='k')
    axs[i][j].axhline(120, color='k')
    axs[i][j].tick_params(labelrotation=0)
    j += 1
    if j%3==0:
        j = 0
        i += 1
    d_methodQuality[method] = np.logical_and(((reconst + ureconst) > 80), ((reconst - ureconst) < 120)).sum()

for x in range(0, 3):
    axs[x][0].set_ylabel('Mass reconstructed [%]')
    axs[-1][x].set_xlabel('Date')
    
plt.show()
print(d_methodQuality)


# %%
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
    fig.savefig(f'correlations/correlation_plots_{key1}.png')
    plt.close()

        

# %%
#Prepare correlation plots whithout outliers



for key1 in ['Cd']:
    i, j = 0, 0
    fig, axs = plt.subplots(6, 8, figsize=(45,30))
    for key2 in list(matrix.keys())[:-7]:
        matrix_nooutliers = matrix[[key1, key2]]
        #matrix_nooutliers['zscore'] = matrix_nooutliers[(np.abs(zscore(matrix_nooutliers)) < 3).all(axis=1)]
        axs[i][j].plot(matrix_nooutliers[key1].where(matrix_nooutliers['Cd']<0.008), matrix_nooutliers[key2], 'o')
        axs[i][j].set_xlabel(key1)
        axs[i][j].set_ylabel(key2)
        j += 1
        if j%8 == 0:
            j=0
            i += 1
        del matrix_nooutliers
    fig.savefig(f'correlations_no_outliers/correlation_plots_{key1}.png')
    plt.show()
    plt.close()

        

# %%
method = 'Simon_2011'
resultNormal = estimation_om_oc(matrix.where(events["Event"]=='no'), method=method)
resultEvent = estimation_om_oc(matrix.where(events["Event"].isin(["S", "SP", "SN"])), method=method)
#print(result.summary())

matrix.where(events['Event'].isin(["S", "SP", "SN"]))["PM2.5"].plot(style='o')

#warm_season_index = matrix.index.where(matrix.index.month >= 9).dropna()
#result = estimation_om_oc(matrix.loc[warm_season_index])
print("No events")
print(resultNormal.summary().as_latex())

print("Events")
print(resultEvent.summary())



# %%
#print(result.fittedvalues)

pred_ols = result.get_prediction()
iv_l = pred_ols.summary_frame()["obs_ci_lower"]
iv_u = pred_ols.summary_frame()["obs_ci_upper"]

matrix_dropped_na = matrix.dropna(axis=0).reset_index(drop=True)
x = matrix_dropped_na.index
y = matrix_dropped_na['PM2.5']

unc_dropped_na = unc.dropna(axis=0).reset_index(drop=True)

with plt.style.context('ggplot'):

    fig, ax = plt.subplots(figsize=(13, 6))

    ax.errorbar(x, y, yerr=unc_dropped_na["PM2.5"],
                color='b', marker='.', capsize=2, capthick=2,
                label="data")
    ax.plot(x, result.fittedvalues, "ro-", label="OLS")
#    ax.plot(x, y, "o", label="data")
    #ax.plot(x, y_true, "b-", label="True")
    #ax.plot(y, result.fittedvalues, "o", label="OLS")
    #lims = [np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    #        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    #        ]
    #ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)

    ax.plot(x, iv_u, "r--", linewidth = 0.75, label='IC')
    ax.plot(x, iv_l, "r--", linewidth = 0.75)
    ax.legend(loc="best")
    ax.set_xlabel('Date')
    ax.set_ylabel('PM2.5 (µg/m$^3$)')
    plt.show()

    fig, ax = plt.subplots(figsize=(13,6))
    ax.plot(x, result.fittedvalues/y * 100, "b.-")
    ax.set_xlabel("Date")
    ax.set_ylabel("Mass explained (%)")
    ax.axhline(80)
    ax.axhline(120)
    plt.show()


#explained_percentage = is_explained.sum()

# %%
fig, ax = plt.subplots()
ax.plot(matrix['SO4'], matrix['Na no sol'].where(~events['Event'].isin(['S', 'SP', 'SN'])), 'o')
ax.plot(matrix['SO4'], matrix['Na no sol'].where(events['Event'].isin(['S', 'SP', 'SN'])), 'd')
ax.set_xlabel('SO$_4$')
ax.set_ylabel('Na total')
plt.show()

fig, ax = plt.subplots()
ax.plot(matrix['SO4'], matrix['Na sol'].where(~events['Event'].isin(['S', 'SP', 'SN'])), 'o')
ax.plot(matrix['SO4'], matrix['Na sol'].where(events['Event'].isin(['S', 'SP', 'SN'])), 'd')
ax.set_xlabel('SO$_4$')
ax.set_ylabel('Na$^+$')
plt.show()

fig, ax = plt.subplots()
ax.plot(matrix['NH4'], matrix['NO3'].where(~events['Event'].isin(['S', 'SP', 'SN'])), 'o')
ax.plot(matrix['NH4'], matrix['NO3'].where(events['Event'].isin(['S', 'SP', 'SN'])), 'd')
ax.set_xlabel('NH$_4$')
ax.set_ylabel('NO$_3$')
plt.show()

fig, ax = plt.subplots()
ax.plot(matrix['PM2.5'], matrix['EC'].where(~events['Event'].isin(['S', 'SP', 'SN'])), 'o')
ax.plot(matrix['PM2.5'], matrix['EC'].where(events['Event'].isin(['S', 'SP', 'SN'])), 'd')
ax.set_xlabel('PM$_{2.5}$')
ax.set_ylabel('EC')
plt.show()

fig, ax = plt.subplots()
ax.plot(matrix['PM2.5'], matrix['OC'].where(~events['Event'].isin(['S', 'SP', 'SN'])), 'o')
ax.plot(matrix['PM2.5'], matrix['OC'].where(events['Event'].isin(['S', 'SP', 'SN'])), 'd')
ax.set_xlabel('PM$_{2.5}$')
ax.set_ylabel('OC')
plt.show()

# %%
## SSA SO4
x = matrix['SO4'].values
y = matrix['Na sol'].values
mask = ~np.isnan(x) & ~np.isnan(y)


slope, intercept, r, p, se  = linregress(x[mask], y[mask])

fig, ax = plt.subplots()
ax.plot(x, y , 'o')
ax.plot(x, intercept + slope * x, label=f'{slope.round(2)} [SO$_{{{4}}}$] + {intercept.round(2)}')
ax.set_xlabel('SO4')
ax.set_ylabel('Na sol')
ax.legend()
fig.savefig('SO4_NAsol.png')
plt.show()

fig, ax = plt.subplots()
ax.plot(matrix['SO4'], matrix['Na no sol'], 'o', label='Na (total)')
ax.plot(x, intercept + slope * x,
       label=f'{slope.round(2)} [SO$_{{{4}}}$] + {intercept.round(2)}', zorder=3)
ax.plot(matrix['SO4'], matrix['Na sol'], 'rx', label='Na$^+$')
ax.set_xlabel('SO4')
ax.set_ylabel('Na total')
ax.legend()
fig.savefig('SO4_NAtotal_withregress.png')
plt.show()

# %%
# SO4 and Cl
fig, ax = plt.subplots()
ax.plot(matrix['SO4'], )

# %%
#Separating K source
#It is often assumed that crustal K is 0.6 Fe

fig, ax = plt.subplots()
ax.plot(matrix['K'], matrix['Na sol'], 'o')
#ax.plot(matrix['K'], matrix['K'] / 0.6)
ax.set_xlabel('K')
ax.set_ylabel('Na$^+$')
plt.show()

fig, ax = plt.subplots()
ax.plot(matrix['K'], matrix['Na sol'].where(matrix['Na sol'] < 0.8), 'o')
ax.set_xlabel('K')
ax.set_ylabel('Na$^+$')
plt.show()

fig, ax = plt.subplots()
ax.plot(matrix['K'], matrix['Na no sol'], 'o')
#ax.plot(matrix['K'], matrix['K'] / 0.6)
ax.set_xlabel('K')
ax.set_ylabel('Na total')
plt.show()

fig, ax = plt.subplots()
ax.plot(matrix['K'], matrix['K'], 'o')
ax.plot(matrix['K'], matrix['Fe'] * 0.6, 'x')
ax.set_xlabel('K')
ax.set_ylabel('0.6 * Fe')
plt.show()

fig, ax = plt.subplots()
ax.plot()

# %%
fig, ax = plt.subplots()
ax.plot(gases['SO4'], gases['SO2'], 'o')
ax.plot(gases['SO4'].where((gases['Na no sol']>2) & (gases['SO4'] <2)), gases['SO2'], 'o')

ax.set_xlabel('SO$_4$')
ax.set_ylabel('SO$_2$')
plt.show()

# %%
display(matrix.loc[matrix['SO4'] > 2])


