# %% plot 3 paneles con funcion
import os
os.chdir('/home/usuario/mdiaz/Documents/paper_facu/')
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
from funciones_pmfBA import plot3panels
from load_data import load_data
from sklearn.metrics import mean_squared_error
import pint

plt.style.use('seaborn-v0_8-paper')
matrix, unc, meteo, gases, events = load_data('data/PMF_BA_fullv4.xlsx', 
                                              'data/PMF_BA_fullv4.xlsx',
                                              'gases_mean.csv', 'data/datos_meteo_blhera5.csv',
                                              'BA_events_testMnew2.xlsx')


# methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994',
#            'Malm_1994', 'Chow_1996', 'Andrews_2000',
#            'Malm_2000', 'Maenhaut_2002', 'DeBell_2006',
#            'Hand_2011','Simon_2011']
methods = ['Simon_2011']

event_columnname="Event_F"
event_labels= ["SI" ,"SF","SO"] #,"SP", "SN" "SL", "S", "SC", "SO"

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
# plot_graph=True
# plot_graphmodall = True
# plot_graph1panel = True
plot_graph=True
plot_graphmodall = False
plot_graph1panel = False
slope={}
stderr={}
intercept={}
intercept_stderr={}
for method in methods:
    # print(method)
    # reconstruccion masica para el metodo original
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

    if (plot_graphmodall==True):
        uncertainty = mass[3]
        total_reconst_mass = mass[0]
        utotal_reconst_mass = mass[2]
        mass = mass[1]
        plot3panels(mass, matrix, uncertainty, unc, total_reconst_mass, utotal_reconst_mass, events,
                event_columnname=event_columnname, event_labels=event_labels, method = method, suffix= 'modall')


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
    #total_reconst_mass, mass, utotal_reconst_mass, uncertainty
    if (plot_graph == True):
         uncertainty = mass[3]
         total_reconst_mass = mass[0]
         utotal_reconst_mass = mass[2]
         mass = mass[1]
         plot3panels(mass, matrix, uncertainty, unc, total_reconst_mass, utotal_reconst_mass, events,
                event_columnname=event_columnname, event_labels=event_labels, method = method, suffix= 'events')

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
# %% Table 3 (chequeado al 18-09-2024)
method = "Simon_2011"
#Estimo beta para no evento
resultNormal = linear_estimation_om_oc(matrix.where(events[event_columnname] == 'no'), method=method, ssa_as_Na=False, display_latex=False)
# Estimo beta para evento
resultEvent = linear_estimation_om_oc(matrix.where(events[event_columnname].isin(event_labels)), method=method,
    ssa_as_Na=False, display_latex=False)

[total_reconst_mass, mass, utotal_reconst_mass, uncertainty] = mass_reconstruction_all(matrix, unc, events, equation=method, 
                                                                                       betas_event=[resultEvent.intercept,resultEvent.slope,1,1,1], 
                                                                                       betas_noevent=[resultNormal.intercept, resultNormal.slope,1,1,1], 
                                                                                       type_reconstruction="events")
mass=pd.DataFrame.from_dict(mass)
mass = mass.rename(columns={'organic_mass': 'Organic mass', 
                            'geological_minerals': 'Geological minerals',
                            'inorganic_ions': 'Inorganic ions', 
                            'elemental_C':'Elemental carbon',
                            'salt': 'Sea salt', 
                            'others': 'KNON'})
keys = ['Organic mass', 'Geological minerals','Inorganic ions','Elemental carbon','Sea salt','Others']
print('Category & Total ($\mu$g/m$^{3}$) & No event ($\mu$g/m$^{3}$) & Event ($\mu$g/m$^{3}$) \\\\ \\hline')
    
for key in keys:
    mass_events=mass[key].where(events[event_columnname].isin(event_labels))
    mass_noevents=mass[key].where(~events[event_columnname].isin(event_labels))
    print(f"{key} & {np.min(mass[key]):.2f} & {np.max(mass[key]):.2f} & {np.mean(mass[key]):.2f} & "
      f"{np.min(mass_noevents):.2f} & {np.max(mass_noevents):.2f} & {np.mean(mass_noevents):.2f} & "
      f"{np.min(mass_events):.2f} & {np.max(mass_events):.2f} &{np.mean(mass_events):.2f} \\\\")


# %%
