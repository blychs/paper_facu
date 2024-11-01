import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
from scipy.stats import linregress, spearmanr, zscore
import statsmodels.api as sm
from funciones_pmfBA import mass_reconstruction, mass_reconstruction_mod, percentage_with_err
from funciones_pmfBA import estimation_om_oc, calculate_seasonal, linear_estimation_om_oc
from funciones_pmfBA import print_stats
from load_data import load_data
#from mass_reconstruction_plot import mass_reconstruction_plot, axvlines


plt.style.use('seaborn-v0_8-paper')
matrix, unc, meteo, gases, events = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                              'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                              'BA_events_testM.xlsx')
# %%
d_methodQuality = {}
methods = ["Simon_2011"]
event_columnname = "Event_M"
event_labels = ["S","SP","SN","SL"]
beta_omoc_event = 2.6
beta_omoc_noevent =2
beta_omoc_all =2.3
plt.style.use('seaborn-v0_8-paper')
fig, axs = plt.subplots(4, 3, figsize=(16, 12))
# methods=['Maenhaut_2002','Hand_2011','Simon_2011']
i, j = 0, 0
for method in methods:
    d_methodQuality[method] = 0
    mass = mass_reconstruction_mod(
        matrix, unc, events=events, equation=method,  event_labels=event_labels, 
        event_column=event_columnname, omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, 
        omoc_all=beta_omoc_all, all_together=False)
    reconst = mass[0]/matrix['PM2.5'] * 100
    ureconst = np.sqrt((1/matrix['PM2.5'] * unc['PM2.5']) **
                       2 + (matrix['PM2.5']/mass[0]/mass[0] * mass[2])**2) * 100
    axs[i][j].errorbar(matrix.index, reconst,  yerr=ureconst, capsize=2,
                       capthick=1, marker='.', ecolor='cornflowerblue', zorder=0)
    axs[i][j].plot(matrix.index, reconst.where(
        events[event_columnname] == 'S'), 'o', label='Smoke', zorder=1)
    axs[i][j].plot(matrix.index, reconst.where(events['Event'] ==
                   'SP'), 'D', label='Smoke previous day', zorder=2)
    axs[i][j].plot(matrix.index, reconst.where(
        events[event_columnname] == 'SN'), 's', label='Smoke next day', zorder=3)
    axs[i][j].plot(matrix.index, reconst.where(
        events[event_columnname] == 'SL'), 'o', label='Smoke local', zorder=4)
    axs[i][j].set_title(method)
    axs[i][j].legend(loc=9)
    axs[i][j].axhline(80, color='k')
    axs[i][j].axhline(120, color='k')
    axs[i][j].tick_params(labelrotation=0)
    j += 1
    if j % 3 == 0:
        j = 0
        i += 1
    d_methodQuality[method] = np.logical_and(
        ((reconst + ureconst) > 80), ((reconst - ureconst) < 120)).sum()

for x in range(0, 3):
    axs[x][0].set_ylabel('Mass reconstructed [%]')
    axs[-1][x].set_xlabel('Date')

plt.show()
print(d_methodQuality)


# %% average mass reconstruction
def average_mass_reconstruction(mass_Hand,mass_Maenhaut,mass_Simon):

    categories = {}

    for key in mass_Hand[1].keys():
        categories[key] = (mass_Simon[1][key] + mass_Hand[1][key] + mass_Maenhaut[1][key])/3
    # Hand no tiene Others y Maenhaut tiene other + trace_elements
    categories['others'] = (mass_Simon[1]['others'] + (mass_Maenhaut[1]['others'] + mass_Maenhaut[1]['trace_elements']))/3
    
    ucategories = {}
    
    for key in mass_Hand[3].keys():
        ucategories[key] = np.linalg.norm([mass_Hand[3][key], mass_Maenhaut[3][key], mass_Simon[3][key]], ord=1, axis=0)/3
    
    ucategories["uothers"] = np.linalg.norm(
        [mass_Simon[3]["uothers"], mass_Maenhaut[3]["uothers"], mass_Maenhaut[3]["utrace_elements"]],
        ord=1, axis=0)/3

    closure = sum(categories.values())
    # print(categories)
    uclosure = np.linalg.norm(ucategories.values(), axis=0)
    print(uclosure)
    return closure, categories, uclosure, ucategories

# %%
def average_mass_reconstruction_MS(mass_Maenhaut,mass_Simon):
    keys = ['organic_mass', 'geological_minerals','inorganic_ions','elemental_C','salt','unexplained']
    ukeys = ['uorganic_mass', 'ugeological_minerals','uinorganic_ions','uelemental_C','usalt','uunexplained']

    categories = {}

    for key in keys:
        categories[key] = (mass_Simon[1][key] + mass_Maenhaut[1][key])/2
    # Hand no tiene Others y Maenhaut tiene other + trace_elements
    categories['others'] = (mass_Simon[1]['others'] + (mass_Maenhaut[1]['others'] + mass_Maenhaut[1]['trace_elements']))/3
    
    ucategories = {}

    for key in ukeys:
        ucategories[key] = np.linalg.norm([ mass_Maenhaut[3][key], mass_Simon[3][key]], axis=0)

    closure = (mass_Simon[0] + mass_Maenhaut[0])/2
    uclosure = np.linalg.norm([mass_Simon[2], mass_Maenhaut[2]], axis=0)
    return closure, categories, uclosure, ucategories
# %% Grafico 3 paneles
def select_events(df, events=events):
    return df.where(events[event_columnname].isin(event_labels))


def select_no_events(df, events=events):
    return df.where(~events[event_columnname].isin(event_labels))

def axvlines(ax=None, xs=[0, 1], ymin=0, ymax=1, **kwargs):
    ax = ax or plt.gca()
    for x in xs:
        ax.axvline(x, ymin=ymin, ymax=ymax, **kwargs)

def axvlines(ax=None, xs=[0, 1], ymin=0, ymax=1, **kwargs):
        ax = ax or plt.gca()
        for x in xs:
            ax.axvline(x, ymin=ymin, ymax=ymax, **kwargs)

methods = ["Simon_2011"]
width = 2.5

for method in methods:

    total_reconst_mass, mass, utotal_reconst_mass, ucategories = mass_reconstruction_mod(
        matrix, unc, events, equation=method, event_labels=event_labels,
        event_column=event_columnname, omoc_event=2.6, omoc_noevent=2,
        omoc_all=2.3, all_together=False, betas_all = betas_all, betas_noevent = betas_noevent,
        betas_event = betas_event
    )

    organic_mass_per = percentage_with_err(
        mass['organic_mass'], matrix['PM2.5'], ucategories['uorganic_mass'], unc['PM2.5'])
    inorganic_ions_per = percentage_with_err(
        mass['inorganic_ions'], matrix['PM2.5'],  ucategories['uinorganic_ions'], unc['PM2.5'])
    geological_minerals_per = percentage_with_err(
        mass['geological_minerals'], matrix['PM2.5'], ucategories['ugeological_minerals'], unc['PM2.5'])
    EC_per = percentage_with_err(
        mass['elemental_C'], matrix['PM2.5'], ucategories['uelemental_C'], unc['PM2.5'])
    plt.style.use('seaborn-v0_8-paper')

    reconst = percentage_with_err(val=total_reconst_mass, uval=utotal_reconst_mass,
                                totalval=matrix['PM2.5'], utotalval=unc['PM2.5'])
 

    fig, ax = plt.subplots(nrows=3, figsize=(7, 7.5), sharex=True, dpi=200)


    fig.suptitle(f'Mass reconstruction - {method}')
    # ax.set_title('Mass reconstructed')
    ax[0].errorbar(matrix.index, matrix['PM2.5'], yerr=unc['PM2.5'],
                color='k', capsize=2, capthick=1, lw=1, marker='.', label='Gravimetric mass', zorder=1)
    ax[0].errorbar(matrix.index, total_reconst_mass, yerr=utotal_reconst_mass, color='red',
                capsize=2, capthick=1, lw=1, marker='.', label='Reconstructed mass', zorder=0)
    ax[0].set_ylabel('PM$_{2.5}$ (µg/m$^3$)')
    ax[0].plot(matrix.index, matrix['PM2.5'].where(events[event_columnname].isin(['S', 'SN', 'SP','SL'])) * 0, 'd',

            color='gray', label='Smoke events', zorder=3)
    # ax[0].plot(matrix.index, matrix['PM2.5'] - total_reconst_mass, '.-')
    # ax[0].plot(matrix.index, events[event_columnname].isin(event_labels), 'X')
    ax[0].legend()
    ax[0].set_ylim(bottom=0, top=60)
    #


    axvlines(ax=ax[0], xs=matrix.index.values, color='silver',
            lw=0.5, linestyle='dotted', zorder=0)
    axvlines(ax=ax[1], xs=matrix.index.values, color='silver',
            lw=0.5, linestyle='dotted', zorder=0)
    axvlines(ax=ax[2], xs=matrix.index.values, color='silver',
            lw=0.5, linestyle='dotted', zorder=0)

    ax[1].bar(matrix.index.values, reconst["perc"].values, width,
               yerr=reconst["uperc"],error_kw={'lw': 1, 'capsize': 2, 'capthick': 1,
                                               'ecolor': 'gray', 'marker': '.'})

    ax[1].axhline(100, linestyle=':', color='k')
    ax[1].axhline(100, linestyle=':', color='k')
    ax[1].axhspan(80, 120, alpha=0.3, color='y')
    ax[1].set_ylabel('Reconstructed mass (%)')
    ax[1].set_xlabel('Date')
    handles, labels = ax[1].get_legend_handles_labels()
    ax[1].legend(reversed(handles), reversed(labels), loc=1)
    ax[1].set_ylim(bottom=0, top=300)
    ax[2].axhline(0, color="gray")
    ax[2].errorbar(matrix.index, total_reconst_mass - matrix['PM2.5'],
                yerr=(unc['PM2.5'] + utotal_reconst_mass), linewidth=0,
                color='tab:blue', capsize=2, capthick=1, elinewidth=1,
                marker='o', label='Gravimetric mass', zorder=3)
    ax[2].errorbar(matrix.index,
                (total_reconst_mass - matrix['PM2.5']).where(
                    events[event_columnname].isin(event_labels)),
                yerr=(unc['PM2.5'] + utotal_reconst_mass), linewidth=0,
                color='tab:red', capsize=2, capthick=1, elinewidth=1,
                marker='o', label='Gravimetric mass', zorder=3)
    ax[2].set_ylim(bottom=-17, top=17)
    ax[2].set_xlabel("date")
    ax[2].set_ylabel("reconstructed - gravimetric ($\mu$g/m$^3$)")
    fig.tight_layout()
    plt.subplots_adjust(hspace=.0)
    plt.subplots_adjust(wspace=.0)
    # fig.savefig(f'images/beta_ne18beta_e23/stacked_bar_daily_percentage_testM_resid_{method}.png')