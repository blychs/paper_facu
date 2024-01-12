#%%
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


plt.style.use('seaborn-v0_8-paper')
matrix, unc, meteo, gases, events = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                              'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                              'BA_events_testM.xlsx')


matrix.describe().to_csv('description_statistics_allM.csv')

methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994',
           'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006',
           'Hand_2011', 'Simon_2011']

event_columnname="Event"
event_labels= ["S", "SP", "SN","SL"]
omoc_noevent=[]
omoc_event=[]
omoc_all=[]

for method in methods:
    resultNormal = estimation_om_oc(matrix.where(
        events[event_columnname] == 'no'), method=method, ssa_as_Na=False)
    omoc_noevent.append(resultNormal.params[1])
    resultEvent = estimation_om_oc(matrix.where(
        events[event_columnname].isin(event_labels)), method=method,
        ssa_as_Na=False)
    omoc_event.append(resultEvent.params[1])
    resultAll = estimation_om_oc(matrix, method=method,
        ssa_as_Na=False)
    omoc_all.append(resultAll.params[1])
    #print(result.summary())

    # matrix.where(events[columnName].isin(["S", "SP", "SN","DS"]))[
    #     "PM2.5"].plot(style='o')

    # warm_season_index = matrix.index.where(matrix.index.month >= 9).dropna()
    # result = estimation_om_oc(matrix.loc[warm_season_index])
    # print(f"{method}")
    #print("No events")
    #print(resultNormal.summary())#.as_latex())
    # print(method, '&', f'{resultNormal.params[1]:.4g}', '&',
    #      f'{resultNormal.bse[1]:.2g}', '&',
    #      f'{resultNormal.pvalues[1]:.3g}', '\\\\')
    #print("Events")
    #print(resultEvent.summary())#.as_latex())
    # print(method, '&', f'{resultEvent.params[1]:.4g}', '&',
    #      f'{resultEvent.bse[1]:.2g}', '&',
    #      f'{resultEvent.pvalues[1]:.3g}', '\\\\')
    # print("All")
    # print(resultAll.summary())#.as_latex())
    # print(method, '&', f'{resultAll.params[1]:.4g}', '&',
    #       f'{resultAll.bse[1]:.2g}', '&',
    #       f'{resultAll.pvalues[1]:.3g}', '\\\\')
    print("No Events & ", method, '&', f'{resultNormal.params[1]:.4g}', '&',
         f'{resultNormal.bse[1]:.2g}', '&',
         f'{resultNormal.pvalues[1]:.3g}', '\\\\')
    print("Events & ",method, '&', f'{resultEvent.params[1]:.4g}', '&',
         f'{resultEvent.bse[1]:.2g}', '&',
         f'{resultEvent.pvalues[1]:.3g}', '\\\\')
    print("All together & ",method, '&', f'{resultAll.params[1]:.4g}', '&',
         f'{resultAll.bse[1]:.2g}', '&',
         f'{resultAll.pvalues[1]:.3g}', '\\\\')
# print(f'{np.mean(omoc_noevent):.2g}', f'{np.mean(omoc_event):.2g}',f'{np.mean(omoc_all):.2g}')
beta_omoc_noevent=np.round(np.mean(omoc_noevent),1)
beta_omoc_event=np.round(np.mean(omoc_event),1)
beta_omoc_all=np.round(np.mean(omoc_all),1)
print(beta_omoc_all,beta_omoc_event,beta_omoc_noevent)

# %%
# #%matplotlib widget
mass_Simon = mass_reconstruction_mod(
    matrix, unc, events=events, equation="Simon_2011",  event_labels=event_labels, event_column=event_columnname, omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, all_together=False)
mass_Hand = mass_reconstruction_mod(
    matrix, unc, events=events, equation="Hand_2011", event_labels=event_labels, event_column=event_columnname, omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, all_together=False)
mass_Maenhaut = mass_reconstruction_mod(
    matrix, unc, events=events, equation="Maenhaut_2002", event_labels=event_labels, event_column=event_columnname, omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, all_together=False)

mass = {}

for key in mass_Hand[1].keys():
    mass[key] = (mass_Simon[1][key] + mass_Hand[1]
                 [key] + mass_Maenhaut[1][key])/3
    
uncertainty = {}

for key in mass_Hand[3].keys():
    uncertainty[key] = np.linalg.norm([mass_Hand[3][key], mass_Maenhaut[3][key],
                                      mass_Simon[3][key]], axis=0)


total_reconst_mass = (mass_Simon[0] + mass_Hand[0] + mass_Maenhaut[0])/3
utotal_reconst_mass = np.linalg.norm(
    [mass_Simon[2], mass_Hand[2], mass_Maenhaut[2]], axis=0)
# print(utotal_reconst_mass)

organic_mass_per = percentage_with_err(
    mass['organic_mass'], matrix['PM2.5'], uncertainty['uorganic_mass'], unc['PM2.5'])
inorganic_ions_per = percentage_with_err(
    mass['inorganic_ions'], matrix['PM2.5'], uncertainty['uinorganic_ions'], unc['PM2.5'])
geological_minerals_per = percentage_with_err(
    mass['geological_minerals'], matrix['PM2.5'], uncertainty['ugeological_minerals'], unc['PM2.5'])
EC_per = percentage_with_err(
    mass['elemental_C'], matrix['PM2.5'], uncertainty['uelemental_C'], unc['PM2.5'])
ssa_per = percentage_with_err(
    mass['salt'], matrix['PM2.5'], uncertainty['usalt'], unc['PM2.5'])
others_per = ((mass_Simon[1]['others'] + mass_Maenhaut[1]['others']) /
              2 + mass_Maenhaut[1]['trace_elements']) / total_reconst_mass * 100
plt.style.use('seaborn-v0_8-paper')

reconst = percentage_with_err(val=total_reconst_mass, uval=utotal_reconst_mass,
                              totalval=matrix['PM2.5'], utotalval=unc['PM2.5'])


# smoke_dates = list(matrix.index.where(events[event_columnname] == 'S').dropna())
# print(smoke_dates)


def select_events(df, events=events):
    return df.where(events[event_columnname].isin(event_labels))


def select_no_events(df, events=events):
    return df.where(~events[event_columnname].isin(event_labels))


width = 2.5

fig, ax = plt.subplots(nrows=2, figsize=(7, 5), sharex=True, dpi=200)


def axvlines(ax=None, xs=[0, 1], ymin=0, ymax=1, **kwargs):
    ax = ax or plt.gca()
    for x in xs:
        ax.axvline(x, ymin=ymin, ymax=ymax, **kwargs)


fig.suptitle('Mass reconstruction')
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
#


def axvlines(ax=None, xs=[0, 1], ymin=0, ymax=1, **kwargs):
    ax = ax or plt.gca()
    for x in xs:
        ax.axvline(x, ymin=ymin, ymax=ymax, **kwargs)


axvlines(ax=ax[0], xs=matrix.index.values, color='silver',
         lw=0.5, linestyle='dotted', zorder=0)
axvlines(ax=ax[1], xs=matrix.index.values, color='silver',
         lw=0.5, linestyle='dotted', zorder=0)

ax[1].bar(matrix.index.values, organic_mass_per['perc'].where(
    matrix['Na sol'].notna()).values, width,  label='OM')
ax[1].bar(matrix.index.values, inorganic_ions_per['perc'].values,
          width,  bottom=organic_mass_per['perc'].values, label='II')
ax[1].bar(matrix.index.values, geological_minerals_per['perc'].values, width,
          bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc']).values, label='GM')
ax[1].bar(matrix.index.values, EC_per['perc'].values, width,
          bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + geological_minerals_per['perc']).values, label='EC')
ax[1].bar(matrix.index.values, ssa_per['perc'].values, width,
          error_kw={'lw': 1, 'capsize': 2, 'capthick': 1,
                    'ecolor': 'gray', 'marker': '.'},
          bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + geological_minerals_per['perc'] + EC_per['perc']).values, label='SSA')
ax[1].bar(matrix.index.values, others_per.values, width, yerr=reconst['uperc'],
          error_kw={'lw': 1, 'capsize': 2, 'capthick': 1,
                    'ecolor': 'gray', 'marker': '.'},
          bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] +
                  geological_minerals_per['perc'] + EC_per['perc'] + ssa_per['perc']).values,
          label='Others')
ax[1].axhline(100, linestyle=':', color='k')
ax[1].axhline(100, linestyle=':', color='k')
ax[1].axhspan(80, 120, alpha=0.2, color='y')
ax[1].set_ylabel('Reconstructed mass (%)')
ax[1].set_xlabel('Date')
handles, labels = ax[1].get_legend_handles_labels()
ax[1].legend(reversed(handles), reversed(labels), loc=1)
fig.tight_layout()
plt.subplots_adjust(hspace=.0)
plt.subplots_adjust(wspace=.0)
fig.savefig('images/stacked_bar_daily_percentage_testM.png')
#plt.show()

# %% Calculo de Tabla 4
mass=pd.DataFrame.from_dict(mass)
keys = ['organic_mass', 'geological_minerals','inorganic_ions','elemental_C','salt','unexplained']
print('Category & Total ($\mu$g/m$^{3}$) & No event ($\mu$g/m$^{3}$) & Event ($\mu$g/m$^{3}$) \\\\ \\hline')
for key in keys:
    mass_events=mass[key].where(events[event_columnname].isin(event_labels))
    mass_noevents=mass[key].where(~events[event_columnname].isin(event_labels))
    print(key,'&',np.round(np.min(mass[key]),decimals=2),'--',np.round(np.max(mass[key]),decimals=2),'(',np.round(np.mean(mass[key]),decimals=2),') &', 
          np.round(np.min(mass_noevents),decimals=2) ,'--', np.round(np.max(mass_noevents),decimals=2), '(', np.round(np.mean(mass_noevents),decimals=2),') &',
          np.round(np.min(mass_events),decimals=2) ,'--', np.round(np.max(mass_events),decimals=2), '(', np.round(np.mean(mass_events),decimals=2),') \\\\')
          
#np.round(np.min(mass[key]),decimals=2),'-',np.round(np.max(mass[key]),decimals=2),'(',np.round(np.mean(mass[key]),decimals=2),') \\\\')
#matrix.where(events['Event'].isin(["S", "SP", "SN"]))["PM2.5"].plot(style='o')
# mass.where(events[event_columnname].isin(event_labels))[key]
# %%
