
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
from mass_reconstruction_plot import mass_reconstruction_plot, axvlines


plt.style.use('seaborn-v0_8-paper')
matrix, unc, meteo, gases, events = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                              'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                              'BA_events_testM.xlsx')


matrix.describe().to_csv('description_statistics_allM.csv')

methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994',
           'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006',
           'Hand_2011','Hand_2011_mod','Simon_2011']

event_columnname="Event_M"
event_labels= ["S", "SP", "SN","SL"]
omoc_noevent=[]
omoc_event=[]
omoc_all=[]

for method in methods:
    resultNormal = estimation_om_oc(matrix.where(
        events[event_columnname] == 'no'), method=method, ssa_as_Na=False, display_latex=True)
    omoc_noevent.append(resultNormal.params[1])
    resultEvent = estimation_om_oc(matrix.where(
        events[event_columnname].isin(event_labels)), method=method,
        ssa_as_Na=False, display_latex=True)
    omoc_event.append(resultEvent.params[1])
    resultAll = estimation_om_oc(matrix, method=method,
        ssa_as_Na=False, display_latex=True)
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
# Hand no tiene Others y Maenhaut tiene other + trace_elements
mass['others'] = (mass_Simon[1]['others'] + (mass_Maenhaut[1]['others'] + mass_Maenhaut[1]['trace_elements']))/3

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
others_per = ((mass_Simon[1]['others'] + (mass_Maenhaut[1]['others'] + mass_Maenhaut[1]['trace_elements']))/3)/ total_reconst_mass * 100
plt.style.use('seaborn-v0_8-paper')

reconst = percentage_with_err(val=total_reconst_mass, uval=utotal_reconst_mass,
                              totalval=matrix['PM2.5'], utotalval=unc['PM2.5'])
def select_events(df, events=events):
    return df.where(events[event_columnname].isin(event_labels))


def select_no_events(df, events=events):
    return df.where(~events[event_columnname].isin(event_labels))


width = 2.5

# mass_reconstruction_plot(matrix,events,total_reconst_mass, utotal_reconst_mass, unc, 
#                          event_labels=event_labels,
#                          savefig=True,imagepath='images/stacked_bar_daily_percentage_testM.png',
#                          showplot=True)
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
mass = mass.rename(columns={'organic_mass': 'Organic mass', 
                            'geological_minerals': 'Geological minerals',
                            'inorganic_ions': 'Inorganic ions', 
                            'elemental_C':'Elemental carbon',
                            'salt': 'Sea salt', 
                            'others': 'Others'})
keys = ['Organic mass', 'Geological minerals','Inorganic ions','Elemental carbon','Sea salt','Others']
print('Category & Total ($\mu$g/m$^{3}$) & No event ($\mu$g/m$^{3}$) & Event ($\mu$g/m$^{3}$) \\\\ \\hline')
for key in keys:
    mass_events=mass[key].where(events[event_columnname].isin(event_labels))
    mass_noevents=mass[key].where(~events[event_columnname].isin(event_labels))
    print(f"{key} & {np.min(mass[key]):.2f} -- {np.max(mass[key]):.2f} ({np.mean(mass[key]):.2f}) & "
      f"{np.min(mass_noevents):.2f} -- {np.max(mass_noevents):.2f} ({np.mean(mass_noevents):.2f}) & "
      f"{np.min(mass_events):.2f} -- {np.max(mass_events):.2f} ({np.mean(mass_events):.2f}) \\\\")

# %% Table 3
# #%matplotlib widget

methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011', 'Hand_2011_mod' , 'Simon_2011']

d_methodQuality = {}
d_methodQuality_modall = {}
d_methodQuality_moddis = {}

for method in methods:
    d_methodQuality[method] = 0
    # mass = mass_reconstruction(matrix, unc, equation=method)
    # reconst = mass[0]/matrix['PM2.5'] * 100
    # ureconst = np.sqrt((1/matrix['PM2.5'] * unc['PM2.5']) **
    #                    2 + (matrix['PM2.5']/mass[0]/mass[0] * mass[2])**2) * 100
    # d_methodQuality[method] = np.logical_and(
    #     ((reconst + ureconst) > 80), ((reconst - ureconst) < 120)).sum()
    
    d_methodQuality_modall[method] = 0
    mass = mass_reconstruction_mod(matrix, unc, events, equation=method, all_together=True)
    reconst = mass[0]/matrix['PM2.5'] * 100
    ureconst = np.sqrt((1/matrix['PM2.5'] * unc['PM2.5']) **
                       2 + (matrix['PM2.5']/mass[0]/mass[0] * mass[2])**2) * 100
    d_methodQuality_modall[method] = np.logical_and(
        ((reconst + ureconst) > 80), ((reconst - ureconst) < 120)).sum()
    
    
    d_methodQuality_moddis[method] = 0
    mass = mass_reconstruction_mod(matrix, unc, events, equation=method, all_together=False)
    reconst = mass[0]/matrix['PM2.5'] * 100
    ureconst = np.sqrt((1/matrix['PM2.5'] * unc['PM2.5']) **
                       2 + (matrix['PM2.5']/mass[0]/mass[0] * mass[2])**2) * 100
    d_methodQuality_moddis[method] = np.logical_and(
        ((reconst + ureconst) > 80), ((reconst - ureconst) < 120)).sum()

#method_quality = pd.concat([pd.DataFrame([d]) for d in [d_methodQuality, d_methodQuality_modall, d_methodQuality_moddis]]).T
#method_quality.set_index("Original","Modified all", "Modified disaggregated")
#print(method_quality)
# print(d_methodQuality)
print(d_methodQuality_modall)
print(d_methodQuality_moddis)

# %%
d_methodQuality = {}

plt.style.use('seaborn-v0_8-paper')
fig, axs = plt.subplots(4, 3, figsize=(16, 12))

i, j = 0, 0
for method in methods:
    d_methodQuality[method] = 0
    mass = mass_reconstruction(matrix, unc, equation=method)
    reconst = mass[0]/matrix['PM2.5'] * 100
    ureconst = np.sqrt((1/matrix['PM2.5'] * unc['PM2.5']) **
                       2 + (matrix['PM2.5']/mass[0]/mass[0] * mass[2])**2) * 100
    axs[i][j].errorbar(matrix.index, reconst,  yerr=ureconst, capsize=2,
                       capthick=1, marker='.', ecolor='cornflowerblue', zorder=0)
    axs[i][j].plot(matrix.index, reconst.where(
        events['Event'] == 'S'), 'o', label='Smoke', zorder=1)
    axs[i][j].plot(matrix.index, reconst.where(events['Event'] ==
                   'SP'), 'D', label='Smoke previous day', zorder=2)
    axs[i][j].plot(matrix.index, reconst.where(
        events['Event'] == 'SN'), 's', label='Smoke next day', zorder=3)
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


# %%
