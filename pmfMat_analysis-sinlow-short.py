
#%% Import libraries and data
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
from funciones_pmfBA import average_mass_reconstruction, axvlines
from load_data import load_data

#from mass_reconstruction_plot import mass_reconstruction_plot, axvlines


plt.style.use('seaborn-v0_8-paper')
matrix, unc, meteo, gases, events = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                              'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                              'BA_events_testM.xlsx')


matrix.describe().to_csv('description_statistics_allM.csv')

methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994',
           'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006',
           'Hand_2011','Simon_2011']


event_columnname="Event_M"
event_labels= ["S", "SP", "SN","SL"]
omoc_noevent=[]
omoc_event=[]
omoc_all=[]

#%%
def select_events(df, events=events):
    return df.where(events[event_columnname].isin(event_labels))


def select_no_events(df, events=events):
    return df.where(~events[event_columnname].isin(event_labels))
#%% Estima parametros omoc 
for method in methods:
    resultNormal = estimation_om_oc(matrix.where(
        events[event_columnname] == 'no'), method=method, ssa_as_Na=False, display_latex=False)
    omoc_noevent.append(resultNormal.params[1])
    resultEvent = estimation_om_oc(matrix.where(
        events[event_columnname].isin(event_labels)), method=method,
        ssa_as_Na=False, display_latex=False)
    omoc_event.append(resultEvent.params[1])
    resultAll = estimation_om_oc(matrix, method=method,
        ssa_as_Na=False, display_latex=False)
    omoc_all.append(resultAll.params[1])
    print("No Events & ", method, '&', f'{resultNormal.params[1]:.4g}', '&',
         f'{resultNormal.bse[1]:.2g}', '&',
         f'{resultNormal.pvalues[1]:.3g}', '\\\\')
    print("Events & ",method, '&', f'{resultEvent.params[1]:.4g}', '&',
         f'{resultEvent.bse[1]:.2g}', '&',
         f'{resultEvent.pvalues[1]:.3g}', '\\\\')
    print("All together & ",method, '&', f'{resultAll.params[1]:.4g}', '&',
         f'{resultAll.bse[1]:.2g}', '&',
         f'{resultAll.pvalues[1]:.3g}', '\\\\')
beta_omoc_noevent=np.round(np.mean(omoc_noevent),1)
beta_omoc_event=np.round(np.mean(omoc_event),1)
beta_omoc_all=np.round(np.mean(omoc_all),1)
print(beta_omoc_all,beta_omoc_event,beta_omoc_noevent)


# %% TypeError: loop of ufunc does not support argument 0 of type dict_values which has no callable conjugate method

beta_omoc_noevent=1.9
beta_omoc_event=2.6
beta_omoc_all=2.3
mass_Simon = mass_reconstruction_mod(
    matrix, unc, events=events, equation="Simon_2011",  event_labels=event_labels, event_column=event_columnname, 
    omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, all_together=False)
perc_reconst = percentage_with_err(mass_Simon[0], matrix["PM2.5"], mass_Simon[2], unc["PM2.5"])
d_methodQuality_Simon = np.logical_and(
   ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)
).sum()
mass_Hand = mass_reconstruction_mod(
    matrix, unc, events=events, equation="Hand_2011", event_labels=event_labels, event_column=event_columnname, 
    omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, all_together=False)
perc_reconst = percentage_with_err(mass_Hand[0], matrix["PM2.5"], mass_Hand[2], unc["PM2.5"])
d_methodQuality_Hand = np.logical_and(
   ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)
).sum()
mass_Maenhaut = mass_reconstruction_mod(
    matrix, unc, events=events, equation="Maenhaut_2002", event_labels=event_labels, event_column=event_columnname, 
    omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, all_together=False)
perc_reconst = percentage_with_err(mass_Maenhaut[0], matrix["PM2.5"], mass_Maenhaut[2], unc["PM2.5"])
d_methodQuality_Maenhaut = np.logical_and(
   ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)).sum()

total_reconst_mass, mass, utotal_reconst_mass, uncertainty = average_mass_reconstruction(mass_Hand, mass_Maenhaut, mass_Simon)
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

width = 2.5

fig, ax = plt.subplots(nrows=3, figsize=(7, 7.5), sharex=True, dpi=200)

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

axvlines(ax=ax[0], xs=matrix.index.values, color='silver',
         lw=0.5, linestyle='dotted', zorder=0)
axvlines(ax=ax[1], xs=matrix.index.values, color='silver',
         lw=0.5, linestyle='dotted', zorder=0)
axvlines(ax=ax[2], xs=matrix.index.values, color='silver',
         lw=0.5, linestyle='dotted', zorder=0)

ax[1].bar(matrix.index.values, organic_mass_per['perc'].where(matrix['Na sol'].notna()).values, 
          width,  label='OM')
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
ax[1].axhspan(80, 120, alpha=0.3, color='y')
ax[1].set_ylabel('reconstructed mass (%)')
# ax[1].legend(ncol=3)
# ax[1].set_xlabel('date')
handles, labels = ax[1].get_legend_handles_labels()
ax[1].legend(reversed(handles), reversed(labels), loc=1,ncol=2)

ax[2].axhline(0, color="gray")
ax[2].errorbar(matrix.index, total_reconst_mass - matrix['PM2.5'],
               yerr=(unc['PM2.5'] + utotal_reconst_mass), linewidth=0,
               color='tab:blue', capsize=2, capthick=1, elinewidth=1,
               marker='o', label='no event', zorder=3)
ax[2].errorbar(matrix.index,
               (total_reconst_mass - matrix['PM2.5']).where(
                   events[event_columnname].isin(event_labels)),
               yerr=(unc['PM2.5'] + utotal_reconst_mass), linewidth=0,
               color='tab:red', capsize=2, capthick=1, elinewidth=1,
               marker='o', label='smoke event', zorder=3)
ax[2].set_ylim(bottom=-20, top=20)
ax[2].set_xlabel("date")
ax[2].set_ylabel("reconstructed - gravimetric ($\mu$g/m$^3$)")
ax[2].legend()
fig.tight_layout()
plt.subplots_adjust(hspace=.0)
plt.subplots_adjust(wspace=.0)
fig.savefig('images/stacked_bar_daily_percentage_testM_resid.png')

perc_reconst = percentage_with_err(val=total_reconst_mass, uval=utotal_reconst_mass,
                              totalval=matrix['PM2.5'], utotalval=unc['PM2.5'])
d_methodQuality_average = np.logical_and(
   ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)
).sum()
perc_reconst = percentage_with_err(mass_Simon[0], matrix["PM2.5"], mass_Simon[2], unc["PM2.5"])
d_methodQuality_Simon = np.logical_and(
   ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)
).sum()
print("Maenhaut:", d_methodQuality_Maenhaut,", Hand: ",d_methodQuality_Hand,
      ", Simon: ",d_methodQuality_Simon," and average: ", d_methodQuality_average)
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

#####################################################################################
# %% Table 3
# #%matplotlib widget
beta_omoc_noevent=1.9       
beta_omoc_event=2.6
beta_omoc_all=2.3
methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011',  'Simon_2011']

d_methodQuality = {}
d_methodQuality_modall = {}
d_methodQuality_moddis = {}

for method in methods:
    d_methodQuality[method] = 0
    mass = mass_reconstruction(matrix, unc, equation=method)
    perc_reconst = percentage_with_err(mass[0], matrix["PM2.5"], mass[2], unc["PM2.5"])
    d_methodQuality[method] = np.logical_and(
   ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)).sum()

    d_methodQuality_modall[method] = 0
    mass = mass_reconstruction_mod(matrix, unc, events, equation=method, 
                                   omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, 
                                   all_together=True)
    perc_reconst = percentage_with_err(mass[0], matrix["PM2.5"], mass[2], unc["PM2.5"])
    d_methodQuality_modall[method] = np.logical_and(
   ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)).sum()

    
    d_methodQuality_moddis[method] = 0
    mass = mass_reconstruction_mod(matrix, unc, events, equation=method, 
                                   omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, 
                                   all_together=False)

    perc_reconst = percentage_with_err(mass[0], matrix["PM2.5"], mass[2], unc["PM2.5"])

    d_methodQuality_moddis[method] = np.logical_and(
    ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)
    ).sum()

method="this_Study"
d_methodQuality[method] = 0
mass_Simon = mass_reconstruction(
    matrix, unc,  equation="Simon_2011")
mass_Hand = mass_reconstruction(
    matrix, unc, equation="Hand_2011")
mass_Maenhaut = mass_reconstruction(
    matrix, unc, equation="Maenhaut_2002")

mass = average_mass_reconstruction(mass_Hand,mass_Maenhaut,mass_Simon)

perc_reconst = percentage_with_err(mass[0], matrix["PM2.5"], mass[2], unc["PM2.5"])
d_methodQuality[method] = np.logical_and(
   ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)).sum()

d_methodQuality_modall[method] = 0

mass_Simon = mass_reconstruction_mod(
    matrix, unc, events=events, equation="Simon_2011",  event_labels=event_labels, event_column=event_columnname, omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, all_together=True)
mass_Hand = mass_reconstruction_mod(
    matrix, unc, events=events, equation="Hand_2011", event_labels=event_labels, event_column=event_columnname, omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, all_together=True)
mass_Maenhaut = mass_reconstruction_mod(
    matrix, unc, events=events, equation="Maenhaut_2002", event_labels=event_labels, event_column=event_columnname, omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, all_together=True)

mass=average_mass_reconstruction(mass_Hand,mass_Maenhaut,mass_Simon)
perc_reconst = percentage_with_err(mass[0], matrix["PM2.5"], mass[2], unc["PM2.5"])
d_methodQuality_modall[method] = np.logical_and(
   ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)
).sum()

d_methodQuality_moddis[method] = 0
mass_Simon = mass_reconstruction_mod(
    matrix, unc, events=events, equation="Simon_2011",  event_labels=event_labels, event_column=event_columnname, omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, all_together=False)
mass_Hand = mass_reconstruction_mod(
    matrix, unc, events=events, equation="Hand_2011", event_labels=event_labels, event_column=event_columnname, omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, all_together=False)
mass_Maenhaut = mass_reconstruction_mod(
    matrix, unc, events=events, equation="Maenhaut_2002", event_labels=event_labels, event_column=event_columnname, omoc_event=beta_omoc_event, omoc_noevent=beta_omoc_noevent, omoc_all=beta_omoc_all, all_together=False)

mass = average_mass_reconstruction(mass_Hand,mass_Maenhaut,mass_Simon)
perc_reconst = percentage_with_err(mass[0], matrix["PM2.5"], mass[2], unc["PM2.5"])

d_methodQuality_moddis[method] = np.logical_and(
   ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)
).sum()

method_quality = pd.concat([pd.DataFrame([d]) for d in [d_methodQuality, d_methodQuality_modall, d_methodQuality_moddis]]).T
method_quality.columns = ["Original","Modified all", "Modified disaggregated"]
print(method_quality)

# %% Grafico 3 paneles

methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011',  'Simon_2011']
# methods = ['Simon_2011']
# methods = ["Simon_2011_mod"]
width = 2.5
onlytotals = False
for method in methods:

    total_reconst_mass, mass, utotal_reconst_mass, ucategories = mass_reconstruction_mod(
        matrix, unc, events, equation=method, event_labels=event_labels,
        event_column=event_columnname, omoc_event=2.3, omoc_noevent=1.7,
        omoc_all=2.3, all_together=False)
    if onlytotals==False:
        organic_mass_per = percentage_with_err(
            mass['organic_mass'], matrix['PM2.5'], ucategories['uorganic_mass'], unc['PM2.5'])
        inorganic_ions_per = percentage_with_err(
            mass['inorganic_ions'], matrix['PM2.5'],  ucategories['uinorganic_ions'], unc['PM2.5'])
        geological_minerals_per = percentage_with_err(
            mass['geological_minerals'], matrix['PM2.5'], ucategories['ugeological_minerals'], unc['PM2.5'])
        EC_per = percentage_with_err(
            mass['elemental_C'], matrix['PM2.5'], ucategories['uelemental_C'], unc['PM2.5'])
        try: 
            ssa_per = percentage_with_err( mass['salt'], matrix['PM2.5'], uncertainty['usalt'], unc['PM2.5'])
        except:
            print (method, "no tiene la categoria ssa")
        try:
            others_per = percentage_with_err( mass['others'], matrix['PM2.5'], uncertainty['uothers'], unc['PM2.5'])
        except:
            print (method, "no tiene la categoria others")    
    
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


    axvlines(ax=ax[0], xs=matrix.index.values, color='silver',
            lw=0.5, linestyle='dotted', zorder=0)
    axvlines(ax=ax[1], xs=matrix.index.values, color='silver',
            lw=0.5, linestyle='dotted', zorder=0)
    axvlines(ax=ax[2], xs=matrix.index.values, color='silver',
            lw=0.5, linestyle='dotted', zorder=0)

    ax[1].bar(matrix.index.values, reconst["perc"].values, width,
               yerr=reconst["uperc"],error_kw={'lw': 1, 'capsize': 2, 'capthick': 1,
                                               'ecolor': 'gray', 'marker': '.'})
    # ax[1].bar(matrix.index.values, inorganic_ions_per['perc'].values,
    #         width,  bottom=organic_mass_per['perc'].values, label='II')
    # ax[1].bar(matrix.index.values, geological_minerals_per['perc'].values, width,
    #         bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc']).values, label='GM')
    # ax[1].bar(matrix.index.values, EC_per['perc'].values, width,
    #         bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + geological_minerals_per['perc']).values, label='EC')
    # ax[1].bar(matrix.index.values, ssa_per['perc'].values, width,
    #         error_kw={'lw': 1, 'capsize': 2, 'capthick': 1,
    #                     'ecolor': 'gray', 'marker': '.'},
    #         bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + geological_minerals_per['perc'] + EC_per['perc']).values, label='SSA')
    # ax[1].bar(matrix.index.values, others_per.values, width, yerr=reconst['uperc'],
    #         error_kw={'lw': 1, 'capsize': 2, 'capthick': 1,
    #                     'ecolor': 'gray', 'marker': '.'},
    #         bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] +
    #                 geological_minerals_per['perc'] + EC_per['perc'] + ssa_per['perc']).values,
    #         label='Others')
    # ax[1].bar(matrix.index.values, )
    ax[1].axhline(100, linestyle=':', color='k')
    ax[1].axhline(100, linestyle=':', color='k')
    ax[1].axhspan(80, 120, alpha=0.3, color='y')
    ax[1].set_ylabel('reconstructed mass (%)')
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
    #ax[2].axhline(-2.5, color="green")
    ax[2].set_xlabel("date")
    ax[2].set_ylabel("reconstructed - gravimetric ($\mu$g/m$^3$)")
    fig.tight_layout()
    plt.subplots_adjust(hspace=.0)
    plt.subplots_adjust(wspace=.0)
    fig.savefig(f'images/beta_ne18beta_e23/stacked_bar_daily_percentage_testM_resid_{method}.png')
    #plt.show()

#%% Solo Simon plot

total_reconst_mass, mass, utotal_reconst_mass, uncertainty = mass_reconstruction_mod(matrix, unc, events, 
                                                                                     equation="Simon_2011", 
                                                                                     omoc_event=beta_omoc_event, 
                                                                                     omoc_noevent=beta_omoc_noevent, 
                                                                                     omoc_all=beta_omoc_all, 
                                                                                     all_together=False)
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

width = 2.5

fig, ax = plt.subplots(nrows=3, figsize=(7, 7.5), sharex=True, dpi=200)

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

axvlines(ax=ax[0], xs=matrix.index.values, color='silver',
         lw=0.5, linestyle='dashed', zorder=0)
axvlines(ax=ax[1], xs=matrix.index.values, color='silver',
         lw=0.5, linestyle='dashed', zorder=0)
axvlines(ax=ax[2], xs=matrix.index.values, color='silver',
         lw=0.5, linestyle='dashed', zorder=0)

ax[1].bar(matrix.index.values, organic_mass_per['perc'].where(matrix['Na sol'].notna()).values, 
          width,  label='OM')
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
ax[1].axhspan(80, 120, alpha=0.3, color='y')
ax[1].set_ylabel('reconstructed mass (%)')
# ax[1].legend(ncol=3)
# ax[1].set_xlabel('date')
handles, labels = ax[1].get_legend_handles_labels()
ax[1].legend(reversed(handles), reversed(labels), loc=1,ncol=2)

ax[2].axhline(0, color="gray")
ax[2].errorbar(matrix.index, total_reconst_mass - matrix['PM2.5'],
               yerr=(unc['PM2.5'] + utotal_reconst_mass), linewidth=0,
               color='tab:blue', capsize=2, capthick=1, elinewidth=1,
               marker='o', label='no event', zorder=3)
ax[2].errorbar(matrix.index,
               (total_reconst_mass - matrix['PM2.5']).where(
                   events[event_columnname].isin(event_labels)),
               yerr=(unc['PM2.5'] + utotal_reconst_mass), linewidth=0,
               color='tab:red', capsize=2, capthick=1, elinewidth=1,
               marker='o', label='smoke event', zorder=3)
ax[2].set_ylim(bottom=-20, top=20)
ax[2].set_xlabel("date")
ax[2].set_ylabel("reconstructed - gravimetric ($\mu$g/m$^3$)")
ax[2].legend()
fig.tight_layout()
plt.subplots_adjust(hspace=.0)
plt.subplots_adjust(wspace=.0)
fig.savefig('images/stacked_bar_daily_percentage_Simon_residuosycategories.png')
