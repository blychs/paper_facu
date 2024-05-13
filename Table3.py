# %% Table 3
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
from scipy.stats import linregress, spearmanr, zscore
import statsmodels.api as sm
from funciones_pmfBA import mass_reconstruction, mass_reconstruction_mod, mass_reconstruction_all
from funciones_pmfBA import percentage_with_err
from funciones_pmfBA import estimation_om_oc, calculate_seasonal, linear_estimation_om_oc
from funciones_pmfBA import average_mass_reconstruction, axvlines
from load_data import load_data
from sklearn.metrics import mean_squared_error

plt.style.use('seaborn-v0_8-paper')
matrix, unc, meteo, gases, events = load_data('data/PMF_BA_fullv3.xlsx', 'data/PMF_BA_fullv3.xlsx',
                                              'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                              'BA_events_testM.xlsx')

# datesdrop=['2019-05-24','2019-05-27','2019-05-30','2019-06-02', '2020-03-01','2020-01-31','2019-08-04','2019-08-07','2019-08-10']
# datesdrop=['2020-03-01','2020-01-31']
# matrix=matrix.drop(datesdrop,axis=0)
# events=events.drop(datesdrop,axis=0)
matrix.describe().to_csv('description_statistics_allM.csv')

methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994',
           'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006',
           'Hand_2011','Simon_2011']
# methods = ['Macias_1981']

event_columnname="Event_M"
event_labels= ["S", "SP", "SN","SL","SC"]
# %% loop 
# #%matplotlib widget

# sin c
beta_omoc_noevent=1.8       
beta_omoc_event=2.3
beta_omoc_all=2.1

# #con c
# beta_omoc_noevent=1.6       
# beta_omoc_event=2.4
# beta_omoc_all=2.0

d_methodQuality = {}
d_methodQuality_modall = {}
d_methodQuality_moddis = {}
d_methodQuality_absolute = {}
d_methodQuality_modall_absolute = {}
d_methodQuality_moddis_absolute = {}
d_methodQualityRMSE = {}
d_methodQuality_modallRMSE = {}
d_methodQuality_moddisRMSE = {}
omoc_noevent=[]
omoc_event=[]
omoc_all=[]
for method in methods:
    # print(method)
    # hago la reconstruccion masica para el metodo original
    d_methodQuality[method] = 0
    d_methodQualityRMSE[method] = 0
    mass = mass_reconstruction_all(matrix, unc, events, equation=method, type_reconstruction="original")
    
    perc_reconst = percentage_with_err(mass[0], matrix["PM2.5"], mass[2], unc["PM2.5"])
    
    d_methodQuality[method] = np.logical_and(
   ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)).sum()
    pm25_for_rmse = matrix['PM2.5'].to_frame()
    pm25_for_rmse["reconstructed"] = mass[0] #total_reconst_mass, mass, utotal_reconst_mass, uncertainty
    pm25_for_rmse = pm25_for_rmse.dropna()
    d_methodQualityRMSE[method] = mean_squared_error(pm25_for_rmse['PM2.5'],pm25_for_rmse['reconstructed'], squared=False)
    d_methodQuality_absolute[method] = np.logical_and(
    ((mass[0]-mass[2])< (matrix["PM2.5"]+unc["PM2.5"])),((mass[0] + mass[2]) > (matrix["PM2.5"]- unc["PM2.5"]))
    ).sum()
    
    # estimo beta para alltogether
    resultAll = linear_estimation_om_oc(matrix, method=method, ssa_as_Na=False, display_latex=False)
    omoc_all.append(resultAll.slope)
    d_methodQuality_modall[method] = 0
    d_methodQuality_modallRMSE[method] = 0
    mass = mass_reconstruction_all(matrix, unc, events, equation=method, 
                                   betas_all=[resultAll.intercept, resultAll.slope, 1,1,1], 
                                   type_reconstruction="alltogether")
    pm25_for_rmse = matrix['PM2.5'].to_frame()
    pm25_for_rmse["reconstructed"] = mass[0] #total_reconst_mass, mass, utotal_reconst_mass, uncertainty
    pm25_for_rmse = pm25_for_rmse.dropna()
    d_methodQuality_modallRMSE[method] = mean_squared_error(pm25_for_rmse['PM2.5'],pm25_for_rmse['reconstructed'], squared=False)
    perc_reconst = percentage_with_err(mass[0], matrix["PM2.5"], mass[2], unc["PM2.5"])
    d_methodQuality_modall[method] = np.logical_and(
   ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)).sum()
    d_methodQuality_modall_absolute[method] = np.logical_and(
    ((mass[0]-mass[2])< (matrix["PM2.5"]+unc["PM2.5"])),((mass[0] + mass[2]) > (matrix["PM2.5"]- unc["PM2.5"]))
    ).sum()

    #Estimo beta para no evento
    resultNormal = linear_estimation_om_oc(matrix.where(events[event_columnname] == 'no'), method=method, ssa_as_Na=False, display_latex=False)
    omoc_noevent.append(resultNormal.slope)
    # Estimo beta para evento
    resultEvent = linear_estimation_om_oc(matrix.where(events[event_columnname].isin(event_labels)), method=method,
        ssa_as_Na=False, display_latex=False)
    omoc_event.append(resultEvent.slope)
    
    # print(method, f'{rms:.02f}')

    d_methodQuality_moddis[method] = 0
    d_methodQuality_moddisRMSE[method] = 0
    mass = mass_reconstruction_all(matrix, unc, events, equation=method, 
                                   betas_event=[resultEvent.intercept,resultEvent.slope,1,1,1], 
                                   betas_noevent=[resultNormal.intercept, resultNormal.slope,1,1,1], 
                                   type_reconstruction="events")

    perc_reconst = percentage_with_err(mass[0], matrix["PM2.5"], mass[2], unc["PM2.5"])

    d_methodQuality_moddis[method] = np.logical_and(
    ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)
    ).sum()
    
    pm25_for_rmse = matrix['PM2.5'].to_frame()
    pm25_for_rmse["reconstructed"] = mass[0] #total_reconst_mass, mass, utotal_reconst_mass, uncertainty
    pm25_for_rmse = pm25_for_rmse.dropna()
    d_methodQuality_moddisRMSE[method] = mean_squared_error(pm25_for_rmse['PM2.5'],pm25_for_rmse['reconstructed'], squared=False)
    # print(method, f'{rms:.02f}')
    d_methodQuality_moddis_absolute[method] = np.logical_and(
    ((mass[0]-mass[2])< (matrix["PM2.5"]+unc["PM2.5"])),((mass[0] + mass[2]) > (matrix["PM2.5"]- unc["PM2.5"]))
    ).sum()
    
    
    # print(method, 'alltogether', resultAll.intercept,resultAll.slope)
    # print(method, 'event',resultEvent.intercept,resultEvent.slope)
    # print(method, 'noevent',resultNormal.intercept, resultNormal.slope)
    
    print(method)
    # Grafico un panel
    
    # mass[0] = total_reconst_mass, 
    # mass[1] = mass, 
    # mass[2] = utotal_reconst_mass, 
    # mass[3] = uncertainty
    
    uncertainty = mass[3]
    total_reconst_mass = mass[0]
    utotal_reconst_mass = mass[2]
    mass = mass[1]
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
    others_per = percentage_with_err(
        mass['others'], matrix['PM2.5'], uncertainty['uothers'], unc['PM2.5']) 
    plt.style.use('seaborn-v0_8-paper')

    reconst = percentage_with_err(val=total_reconst_mass, uval=utotal_reconst_mass,
                                totalval=matrix['PM2.5'], utotalval=unc['PM2.5'])

    width = 2.5

    fig, ax = plt.subplots(nrows=1, figsize=(7, 4), sharex=True, dpi=200)
    fig.suptitle(method)
    ax.errorbar(matrix.index, matrix['PM2.5'], yerr=unc['PM2.5'],
                color='k', capsize=2, capthick=1, lw=0.5, linestyle='--', marker='.', 
                label='Gravimetric mass', zorder=4)
    ax.set_ylabel('PM$_{2.5}$ (µg/m$^3$)')
    # ax.plot(matrix.index, matrix['PM2.5'].where(events[event_columnname].isin(event_labels))*0+3, 'd', 
    #            color='gray', label='Smoke events', zorder=8)
    # ax.errorbar(matrix.index, matrix['PM2.5'].where(events[event_columnname].isin(event_labels)), yerr=unc['PM2.5'],
    #                color='r', capsize=2, capthick=1, lw=0.5, linestyle='None', marker='.', label='Smoke events', zorder=5)
    ax.plot(matrix.index, matrix['PM2.5'].where(events[event_columnname].isin(event_labels)), '.',
            color='r', label='Smoke events', zorder=8)

    axvlines(ax=ax, xs=matrix.index.values, color='silver',
            lw=0.5, linestyle='dashed', zorder=0)
    # axvlines(ax=ax, xs=matrix.index.values, color='silver',
    #          lw=0.5, linestyle='dashed', zorder=0)
    # axvlines(ax=ax, xs=matrix.index.values, color='silver',
    #          lw=0.5, linestyle='dashed', zorder=0)

    ax.bar(matrix.index.values, mass['organic_mass'].where(matrix['Na sol'].notna()).values, 
            width,  label='OM')
    ax.bar(matrix.index.values, mass['inorganic_ions'].values,
            width,  bottom=mass['organic_mass'].values, label='II')
    ax.bar(matrix.index.values, mass['geological_minerals'].values, width,
            bottom=(mass['inorganic_ions'] + mass['organic_mass']).values, label='GM')
    ax.bar(matrix.index.values, mass['elemental_C'].values, width,
            bottom=(mass['inorganic_ions'] + mass['organic_mass'] + mass['geological_minerals']).values, 
            label='EC')
    ax.bar(matrix.index.values, mass['salt'].values, width,
            bottom=(mass['elemental_C']+mass['inorganic_ions'] + mass['organic_mass'] + mass['geological_minerals']).values, 
            label='SSA')
    ax.bar(matrix.index.values, mass['others'].values, width, yerr=utotal_reconst_mass,
            error_kw={'lw': 1, 'capsize': 2, 'capthick': 1,
                        'ecolor': '#696462', 'marker': 'o'},
            bottom=(mass['salt']+ mass['elemental_C']+mass['inorganic_ions'] + mass['organic_mass'] + mass['geological_minerals']).values,
            label='Others')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'images/stacked_bar_absolute_{method}.png')
    
    # grafico 3 paneles
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
    # others_per = ((mass_Simon[1]['others'] + (mass_Maenhaut[1]['others'] + mass_Maenhaut[1]['trace_elements']))/3)/ total_reconst_mass * 100
    others_per = (mass['others'])/ total_reconst_mass * 100

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
    ax[0].plot(matrix.index, matrix['PM2.5'].where(events[event_columnname].isin(event_labels)) * 0, 'd',

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
            bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + geological_minerals_per['perc']).values, 
            label='EC')
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
    fig.savefig(f'images/stacked_bar_3panels_{method}.png')
# method="Simon_2011_linmod"

# total_reconst_mass, mass, utotal_reconst_mass, uncertainty = mass_reconstruction_mod(matrix, unc, events, 
#                                                                                      equation="Simon_2011_linmod", 
#                                                                                      omoc_event=beta_omoc_event, 
#                                                                                      omoc_noevent=beta_omoc_noevent, 
#                                                                                      omoc_all=beta_omoc_all, 
#                                                                                      all_together=False,
#                                                                                      betas_event=[resultEvent.slope, resultEvent.intercept,1,1],
#                                                                                      betas_noevent=[resultNormal.slope, resultNormal.intercept,1,1])
# organic_mass_per = percentage_with_err(
#     mass['organic_mass'], matrix['PM2.5'], uncertainty['uorganic_mass'], unc['PM2.5'])
# inorganic_ions_per = percentage_with_err(
#     mass['inorganic_ions'], matrix['PM2.5'], uncertainty['uinorganic_ions'], unc['PM2.5'])
# geological_minerals_per = percentage_with_err(
#     mass['geological_minerals'], matrix['PM2.5'], uncertainty['ugeological_minerals'], unc['PM2.5'])
# EC_per = percentage_with_err(
#     mass['elemental_C'], matrix['PM2.5'], uncertainty['uelemental_C'], unc['PM2.5'])
# ssa_per = percentage_with_err(
#     mass['salt'], matrix['PM2.5'], uncertainty['usalt'], unc['PM2.5'])
# # others_per = ((mass_Simon[1]['others'] + (mass_Maenhaut[1]['others'] + mass_Maenhaut[1]['trace_elements']))/3)/ total_reconst_mass * 100
# others_per = (mass['others'])/ total_reconst_mass * 100

# perc_reconst = percentage_with_err(total_reconst_mass, matrix["PM2.5"], utotal_reconst_mass, unc["PM2.5"])

# d_methodQuality_moddis[method] = np.logical_and(
#    ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)
# ).sum()

# pm25_for_rmse = matrix['PM2.5'].to_frame()
# pm25_for_rmse["reconstructed"] = total_reconst_mass #total_reconst_mass, mass, utotal_reconst_mass, uncertainty
# pm25_for_rmse = pm25_for_rmse.dropna()
# d_methodQuality_moddisRMSE[method] = mean_squared_error(pm25_for_rmse['PM2.5'],pm25_for_rmse['reconstructed'], squared=False)

method_quality = pd.concat([pd.DataFrame([d]) for d in [d_methodQuality, d_methodQuality_modall, d_methodQuality_moddis]]).T
method_quality.columns = ["Original","Modified all", "Modified disaggregated"]
print(method_quality)

method_RMSE = pd.concat([pd.DataFrame([d]) for d in [d_methodQualityRMSE, d_methodQuality_modallRMSE, d_methodQuality_moddisRMSE]]).T
method_RMSE.columns = ["Original","Modified all", "Modified disaggregated"]
print(method_RMSE)

method_absolute = pd.concat([pd.DataFrame([d]) for d in [d_methodQuality_absolute, d_methodQuality_modall_absolute, d_methodQuality_moddis_absolute]]).T
method_absolute.columns = ["Original","Modified all", "Modified disaggregated"]
print(method_absolute)
# # %%

# %%
