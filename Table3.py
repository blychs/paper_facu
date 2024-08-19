# %% Table 3
import os
os.chdir('/home/mdiaz/Documents/paper_facu/')
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
import pint

plt.style.use('seaborn-v0_8-paper')
matrix, unc, meteo, gases, events = load_data('data/PMF_BA_fullv3.xlsx', 
                                              'data/PMF_BA_fullv3.xlsx',
                                              'gases_mean.csv', 'data/datos_meteo_blhera5.csv',
                                              'BA_events_testMnew.xlsx')
#%%
# datesdrop=['2019-05-24','2019-05-27','2019-05-30','2019-06-02', '2020-03-01','2020-01-31','2019-08-04','2019-08-07','2019-08-10']
# datesdrop=['2020-03-01','2020-01-31']
# datesdrop = ['2019-08-16'] # por que habiamos sacado este dia?
# matrix=matrix.drop(datesdrop,axis=0)
# events=events.drop(datesdrop,axis=0)
# unc=unc.drop(datesdrop,axis=0)

matrix.describe().to_csv('description_statistics_allM.csv')

# methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994',
#            'Malm_1994', 'Chow_1996', 'Andrews_2000',
#            'Malm_2000', 'Maenhaut_2002', 'DeBell_2006',
#            'Hand_2011','Simon_2011']
methods = ['Simon_2011','Lichtig_2024']

event_columnname="Event_F"
event_labels= ["SI" ,"SF","SO"] # "SL", "S", "SC", "SO"
# event_columnname="Event_M"
# event_labels= ["S", "SL" ,"SP", "SN","SC"] 
# #%matplotlib widget

# # sin c
# beta_omoc_noevent=1.8       
# beta_omoc_event=2.3
# beta_omoc_all=2.1

# # sin c
# beta_omoc_noevent=1.8       
# beta_omoc_event=2.5
# beta_omoc_all=2.2
# #con c
# beta_omoc_noevent=1.6       
# beta_omoc_event=2.4
# beta_omoc_all=2.0

d_methodQuality = {}
d_methodQuality_modall = {}
d_methodQuality_moddis = {}
d_methodQuality_absolute ={}
d_methodQuality_modall_absolute = {}
d_methodQuality_moddis_absolute = {}
d_methodQualityRMSE = {}
d_methodQuality_modallRMSE = {}
d_methodQuality_moddisRMSE = {}
omoc_noevent=[]
omoc_event=[]
omoc_all=[]
plot_graph=True
plot_graph1panel=False
slope={}
stderr={}
intercept={}
intercept_stderr={}
for method in methods:
    # print(method)
    # hago la reconstruccion masica para el metodo original
    d_methodQuality[method] = 0
    d_methodQualityRMSE[method] = 0
    mass = mass_reconstruction_all(matrix, unc, events, equation=method, type_reconstruction="original")
    df=pd.DataFrame.from_dict(mass[1])
    df.to_csv(method+'_original_mass.csv')
    
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
    slope[method,"all"]=resultAll.slope
    stderr[method,"all"]=resultAll.stderr
    intercept[method,"all"]=resultAll.intercept
    intercept_stderr[method,"all"]=resultAll.intercept_stderr
    d_methodQuality_modall[method] = 0
    d_methodQuality_modallRMSE[method] = 0
    mass = mass_reconstruction_all(matrix, unc, events, equation=method, 
                                   betas_all=[resultAll.intercept, resultAll.slope, 1,1,1], 
                                   type_reconstruction="alltogether")
    df=pd.DataFrame.from_dict(mass[1])
    df.to_csv(method+'_alltogether_mass.csv')
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

    if (plot_graph==True):
        # grafico 3 paneles
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
        # others_per = ((mass_Simon[1]['others'] + (mass_Maenhaut[1]['others'] + mass_Maenhaut[1]['trace_elements']))/3)/ total_reconst_mass * 100
        others_per = (mass['others']+mass['residual']+mass['trace_elements'])/ total_reconst_mass * 100

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

        ax[1].bar(matrix.index.values, organic_mass_per['perc'].where(matrix['Na sol'].notna()).values, width,  
                  label='OM')
        ax[1].bar(matrix.index.values, inorganic_ions_per['perc'].values,
                width,  bottom=organic_mass_per['perc'].values, 
                label='II')
        ax[1].bar(matrix.index.values, geological_minerals_per['perc'].values, width,
                bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc']).values, 
                label='GM')
        ax[1].bar(matrix.index.values, EC_per['perc'].values, width,
                bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + geological_minerals_per['perc']).values, 
                label='EC')
        ax[1].bar(matrix.index.values, ssa_per['perc'].values, width,
                bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + geological_minerals_per['perc'] + EC_per['perc']).values, 
                label='SSA')
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
        ax[1].legend(reversed(handles), reversed(labels), loc=9,ncol=6)

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
        fig.savefig(f'images/stacked_bar_3panels_modall_{method}.png')


    #Estimo beta para no evento
    resultNormal = linear_estimation_om_oc(matrix.where(events[event_columnname] == 'no'), method=method, ssa_as_Na=False, display_latex=False)
    omoc_noevent.append(resultNormal.slope)
    # Estimo beta para evento
    resultEvent = linear_estimation_om_oc(matrix.where(events[event_columnname].isin(event_labels)), method=method,
        ssa_as_Na=False, display_latex=False)
    omoc_event.append(resultEvent.slope)

    stderr[method, "event"]=resultEvent.stderr
    slope[method,"no event"]=resultNormal.slope
    stderr[method,"no event"]=resultNormal.stderr
    intercept[method,"no event"]=resultNormal.intercept
    intercept_stderr[method,"no event"]=resultNormal.intercept_stderr
    slope[method,"event"]=resultEvent.slope
    stderr[method,"event"]=resultEvent.stderr
    intercept[method,"event"]=resultEvent.intercept
    intercept_stderr[method,"event"]=resultEvent.intercept_stderr    
    # print(method, f'{rms:.02f}')

    d_methodQuality_moddis[method] = 0
    d_methodQuality_moddisRMSE[method] = 0
    mass = mass_reconstruction_all(matrix, unc, events, equation=method, 
                                   betas_event=[resultEvent.intercept,resultEvent.slope,1,1,1], 
                                   betas_noevent=[resultNormal.intercept, resultNormal.slope,1,1,1], 
                                   type_reconstruction="events")
    df=pd.DataFrame.from_dict(mass[1])
    df.to_csv(method+'_events_mass.csv')

    perc_reconst = percentage_with_err(mass[0], matrix["PM2.5"], mass[2], unc["PM2.5"])

    d_methodQuality_moddis[method] = np.logical_and(
    ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)
    ).sum()
    prueba=np.logical_and(
    ((perc_reconst["perc"] + perc_reconst["uperc"]) > 80),((perc_reconst["perc"] - perc_reconst["uperc"]) < 120)
    )
    prueba2=perc_reconst["perc"] + perc_reconst["uperc"]
    prueba3=perc_reconst["perc"] - perc_reconst["uperc"]
#     print(prueba)
#     print(prueba2)
#     print(prueba3)
    pm25_for_rmse = matrix['PM2.5'].to_frame()
    pm25_for_rmse["reconstructed"] = mass[0] #total_reconst_mass, mass, utotal_reconst_mass, uncertainty
    pm25_for_rmse = pm25_for_rmse.dropna()
    d_methodQuality_moddisRMSE[method] = mean_squared_error(pm25_for_rmse['PM2.5'],pm25_for_rmse['reconstructed'], squared=False)
    # print(method, f'{rms:.02f}')
    d_methodQuality_moddis_absolute[method] = np.logical_and(
    ((mass[0]-mass[2])< (matrix["PM2.5"]+unc["PM2.5"])),((mass[0] + mass[2]) > (matrix["PM2.5"]-unc["PM2.5"]))
    ).sum()
    
    
    print(method, 'alltogether', resultAll.intercept,resultAll.slope)
    print(method, 'event',resultEvent.intercept,resultEvent.slope)
    print(method, 'noevent',resultNormal.intercept, resultNormal.slope)
    
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
    if (plot_graph1panel == True):
        fig, ax = plt.subplots(nrows=1, figsize=(7, 4), sharex=True, dpi=200)
        fig.suptitle(method)
        ax.errorbar(matrix.index, matrix['PM2.5'], yerr=unc['PM2.5'],
                    color='k', capsize=2, capthick=1, lw=0.5, linestyle='--', marker='.', 
                    label='Gravimetric mass', zorder=4)
        ax.set_ylabel('PM$_{2.5}$ (µg/m$^3$)')
        ax.plot(matrix.index, matrix['PM2.5'].where(events[event_columnname].isin(event_labels)), '.',
                color='r', label='Smoke events', zorder=8)

        axvlines(ax=ax, xs=matrix.index.values, color='silver',
                lw=0.5, linestyle='dashed', zorder=0)
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
#     others_per = (mass['others']+mass['residual']+mass['trace_elements'])/ total_reconst_mass * 100
    others_per = (mass['others']+mass['residual']+mass['trace_elements'])/ matrix['PM2.5'] * 100

    plt.style.use('seaborn-v0_8-paper')

    reconst = percentage_with_err(val=total_reconst_mass, uval=utotal_reconst_mass,
                                totalval=matrix['PM2.5'], utotalval=unc['PM2.5'])

    width = 2.5
    if (plot_graph == True):
        fig, ax = plt.subplots(nrows=3, figsize=(7, 7.5), sharex=True, dpi=200)

        ax[0].errorbar(matrix.index, matrix['PM2.5'], yerr=unc['PM2.5'],
                    color='k', capsize=2, capthick=1, lw=1, marker='.', label='gravimetric mass', zorder=1)
        ax[0].errorbar(matrix.index, total_reconst_mass, yerr=utotal_reconst_mass, color='red',
                    capsize=2, capthick=1, lw=1, marker='.', label='reconstructed mass', zorder=0)
        ax[0].set_ylabel('PM$_{2.5}$ (µg/m$^3$)')
        ax[0].set_ylim(bottom=0, top=62)
        ax[0].plot(matrix.index, matrix['PM2.5'].where(events[event_columnname].isin(event_labels)) * 0+2, 'd',

                color='gray', label='smoke event', zorder=5)
        handles, labels = ax[0].get_legend_handles_labels()
        ax[0].legend(reversed(handles), reversed(labels),loc=9,ncol=3)

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
                  bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] + 
                          geological_minerals_per['perc'] + EC_per['perc']).values, 
                  label='SSA')
        ax[1].bar(matrix.index.values, others_per.values, width, yerr=reconst['uperc'],
                error_kw={'lw': 1, 'capsize': 2, 'capthick': 1,
                            'ecolor': 'gray', 'marker': '.'},
                bottom=(inorganic_ions_per['perc'] + organic_mass_per['perc'] +
                        geological_minerals_per['perc'] + EC_per['perc'] + ssa_per['perc']).values,
                label='others')
        ax[1].axhline(100, linestyle=':', color='k')
        ax[1].axhline(100, linestyle=':', color='k')
        ax[1].axhspan(80, 120, alpha=0.3, color='y')
        ax[1].set_ylabel('reconstructed mass (%)')
        ax[1].set_ylim(bottom=0, top=340)
        handles, labels = ax[1].get_legend_handles_labels()
        ax[1].legend(loc=9,ncol=6)

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
        ax[2].set_ylim(bottom=-15, top=17)
        ax[2].set_xlabel("date")
        ax[2].set_ylabel("residuals ($\mu$g/m$^3$)")
        ax[2].legend(loc=9,ncol=2)
        ax[0].text(0.01, 0.95, '(a)', transform=ax[0].transAxes, fontsize=12, verticalalignment='top')
        ax[1].text(0.01, 0.95, '(b)', transform=ax[1].transAxes, fontsize=12, verticalalignment='top')
        ax[2].text(0.01, 0.95, '(c)', transform=ax[2].transAxes, fontsize=12, verticalalignment='top')

        fig.tight_layout()
        plt.subplots_adjust(hspace=.0)
        plt.subplots_adjust(wspace=.0)
        fig.savefig(f'images/stacked_bar_3panels_{method}.png')

# plt.close('all')

method_quality = pd.concat([pd.DataFrame([d]) for d in [d_methodQuality, d_methodQuality_modall, d_methodQuality_moddis]]).T
method_quality.columns = ["Original","Modified all", "Modified disaggregated"]
# print(method_quality)

method_RMSE = pd.concat([pd.DataFrame([d]) for d in [d_methodQualityRMSE, d_methodQuality_modallRMSE, d_methodQuality_moddisRMSE]]).T
method_RMSE.columns = ["Original","Modified all", "Modified disaggregated"]
method_RMSE = method_RMSE.round(2)
# print(method_RMSE)

method_absolute = pd.concat([pd.DataFrame([d]) for d in [d_methodQuality_absolute, d_methodQuality_modall_absolute, d_methodQuality_moddis_absolute]]).T
method_absolute.columns = ["Original","Modified all", "Modified disaggregated"]
# print(method_absolute)

method_quality = pd.concat([method_RMSE, method_quality], axis=1)

latex_table = method_quality.to_latex(escape=False, index_names=False, float_format="%.2f")

# Reemplazar \toprule, \middlerule y \bottomrule por \hline
latex_table = latex_table.replace('\\toprule', '\\hline')
latex_table = latex_table.replace('\\midrule', '\\hline')
latex_table = latex_table.replace('\\bottomrule', '\\hline')

# Imprimir la tabla en formato LaTeX
print(latex_table)
# print(method_quality)
# print method method_quality y RMSE

table1 = pd.concat([pd.DataFrame([d]) for d in [slope, stderr, intercept, 
                                                intercept_stderr]]).T
table1.columns = ["slope","stderr", "intercept", "intercept_stderr"]
# index = pd.MultiIndex.from_tuples(data.keys(), names=['source', 'type'])
# table1 = pd.DataFrame(list(data.values()), index=index, columns=["slope", "stderr", "intercept", "intercept_stderr"])
table1.index.names = ['source', 'type']
table1 = table1.unstack(level='type')

# Reorganizar las columnas
ordered_columns = [
    ('slope', 'no event'), ('stderr', 'no event'), ('intercept', 'no event'), ('intercept_stderr', 'no event'),
    ('slope', 'event'), ('stderr', 'event'), ('intercept', 'event'), ('intercept_stderr', 'event'),
    ('slope', 'all'), ('stderr', 'all'), ('intercept', 'all'), ('intercept_stderr', 'all')
]
table1 = table1[ordered_columns]

# Aplanar el MultiIndex en las columnas
table1.columns = ['_'.join(col) for col in table1.columns]

# Exportar a LaTeX
latex_table = table1.to_latex(escape=False, index=True, float_format="%.2f")

print(latex_table)
latex_table2 = table1.to_latex(escape=False, index_names=False, float_format="%.2f")



#%%
class Compuesto:
    def __init__(self, nombre, columna, tipo, unidad_entrada, unidad_salida,formato, category, decimals):
        self.nombre = nombre
        self.columna = columna
        self.tipo = tipo
        self.unidad_entrada = unidad_entrada
        self.unidad_salida = unidad_salida
        self.formato = formato
        self.category = category
        self.decimals = decimals
class DataFrameToLatex:
    def __init__(self, data, compuestos):
        """
        Inicializa la instancia con un DataFrame y una lista de compuestos.
        
        :param data: Diccionario de datos para crear el DataFrame.
        :param compuestos: Lista de objetos Compuesto.
        """
        self.df = pd.DataFrame(data)
        self.compuestos = compuestos

    def format_columns(self):
        """
        Aplica diferentes formatos a diferentes columnas según la configuración de los compuestos.
        """
        # for compuesto in self.compuestos:
        #     # Aplica el formato especificado en la configuración del compuesto
        #     if compuesto.formato == "float2":
        #         self.df[compuesto.nombre] = self.df[compuesto.nombre].map(lambda x: f"{x:.2f}")
        #     elif compuesto.formato == "float3":
        #         self.df[compuesto.nombre] = self.df[compuesto.nombre].map(lambda x: f"{x:.3f}")
        #     elif compuesto.formato == "float1":
        #         self.df[compuesto.nombre] = self.df[compuesto.nombre].map(lambda x: f"{x:.1f}")
        #     elif compuesto.formato == "integer":
        #         self.df[compuesto.nombre] = self.df[compuesto.nombre].map(lambda x: f"{x:d}")

    def format_columns(self):
        """
        Aplica diferentes formatos a diferentes columnas según la configuración de los compuestos.
        """
        for compuesto in self.compuestos:
            col_name = compuesto.nombre
            if col_name in self.df.columns:
                # Asegurarse de que la columna sea numérica antes de aplicar el formato
                if pd.api.types.is_numeric_dtype(self.df[col_name]):
                    # Aplica el formato especificado en la configuración del compuesto
                    if compuesto.formato == "float2":
                        self.df[col_name] = self.df[col_name].map(lambda x: f"{x:.2f}")
                    elif compuesto.formato == "float3":
                        self.df[col_name] = self.df[col_name].map(lambda x: f"{x:.3f}")
                    elif compuesto.formato == "float1":
                        self.df[col_name] = self.df[col_name].map(lambda x: f"{x:.1f}")
                    elif compuesto.formato == "integer":
                        self.df[col_name] = self.df[col_name].map(lambda x: f"{x:.0f}")
                else:
                    print(f"Warning: Column {col_name} is not numeric and will not be formatted.")

    def to_latex(self, escape=False, index_names=False):
        """
        Convierte el DataFrame formateado a una tabla LaTeX.

        :param escape: Booleano para escapar caracteres LaTeX especiales.
        :param index_names: Booleano para incluir/excluir los nombres de los índices.
        :return: Cadena con la tabla LaTeX.
        """
        return self.df.to_latex(escape=escape, index_names=index_names)
    
# Cargar la información de los compuestos desde el archivo CSV
info_compuestos_df = pd.read_csv('data/info_compuestos.csv')
compuestos = []
for _, row in info_compuestos_df.iterrows():
    compuesto = Compuesto(row['nombre'], row['columna'], row['tipo'], row['unidad_entrada'], row['unidad_salida'], row['formato'], row['category'], row['decimals'])
    compuestos.append(compuesto)
# Crear un registro de unidades
ureg = pint.UnitRegistry()

# Suponiendo que tienes un DataFrame llamado matrix y una lista de objetos Compuesto llamada compuestos

# Reemplazar los valores en el DataFrame con las unidades específicas
for compuesto in compuestos:
    # print(compuesto.nombre)
    # Verificar si la columna del compuesto está presente en el DataFrame
    if compuesto.nombre in matrix.columns:
        # Obtener la función de conversión de unidades
        conversion_factor= ureg(compuesto.unidad_entrada).to(compuesto.unidad_salida)
        # Aplicar la función de conversión a la columna correspondiente del DataFrame
        matrix[compuesto.nombre] = matrix[compuesto.nombre]*conversion_factor if pd.notna(
            matrix[compuesto.nombre]).all() else matrix[compuesto.nombre]

# %%        
matrix_seasonal = calculate_seasonal(matrix)
matrix_seasonal = matrix_seasonal
matrix_seasonal = matrix_seasonal.drop([
    "ECPk1 C", "ECPk2 C", "ECPk3 C", "ECPk4 C",
    "ECPk5 C", "ECPk6 C", "Na total", "OCPk1 C",
    "OCPk2 C", "OCPk3 C", "OCPk4 C", "OCPk5 C",
    "Pyrol C", "Na no sol"], axis=0)

matrix_seasonal['Tipo'] = matrix_seasonal.index.map(lambda x: 
        next((compuesto.tipo for compuesto in compuestos if compuesto.nombre == x), None))

# # Verificar los tipos mapeados en matrix_seasonal
# print("Tipos mapeados en matrix_seasonal:")
# print(matrix_seasonal['Tipo'])


# Ordenar el DataFrame por el tipo
matrix_seasonal_reset = matrix_seasonal.reset_index()
matrix_seasonal_sorted = matrix_seasonal_reset.sort_values(by=['Tipo', 'index'])
matrix_seasonal_sorted = matrix_seasonal_sorted.set_index('index')
matrix_seasonal_sorted.index.name = matrix_seasonal.index.name
matrix_seasonal = matrix_seasonal_sorted
matrix_seasonal = matrix_seasonal.drop(columns=['Tipo'])
mapeo_nombre_a_columna = {compuesto.nombre: compuesto.columna for compuesto in compuestos}
matrix_seasonal = matrix_seasonal.rename(index=mapeo_nombre_a_columna)

# Pm2.5 round2
# Al ng
# Cu  
#%%
# Crear una instancia de DataFrameToLatex y formatear las columnas
df_to_latex = DataFrameToLatex(matrix_seasonal, compuestos)
df_to_latex.format_columns()

# Convertir el DataFrame a una tabla LaTeX
latex_table = df_to_latex.to_latex(escape=False, index_names=False)

print(latex_table)
#%%
# for compuesto in compuestos:
#     matrix_seasonal[compuesto.columna] = matrix_seasonal[compuesto.columna].map(lambda x: f"{x:.2f}")

# latex_table = matrix_seasonal.to_latex(escape=False, index_names=False)

latex_table = matrix_seasonal.to_latex(escape=False, index_names=False, float_format="%.3g")

# Reemplazar \toprule, \middlerule y \bottomrule por \hline
latex_table = latex_table.replace('\\toprule', '\\hline')
latex_table = latex_table.replace('\\midrule', '\\hline')
latex_table = latex_table.replace('\\bottomrule', '\\hline')

# Imprimir la tabla en formato LaTeX
print(latex_table)
# %%
pm_neg = (matrix["PM2.5"] - mass[1]["inorganic_ions"]
          - mass[1]["salt"] 
          - mass[1]["geological_minerals"]
          - mass[1]["elemental_C"]
          - mass[1]["others"]) 
fig, ax = plt.subplots()
ax.plot(matrix["OC"], pm_neg, "ob")
ax.plot(matrix["OC"].where(events[event_columnname].isin(event_labels)),
        pm_neg.where(events[event_columnname].isin(event_labels)), "or")
ax.set_xlabel("OC")
ax.set_ylabel("PM2.5 - cat")
plt.show()
# %%
