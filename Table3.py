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
d_methodQualityRMSE = {}
d_methodQuality_modallRMSE = {}
d_methodQuality_moddisRMSE = {}
omoc_noevent=[]
omoc_event=[]
omoc_all=[]
for method in methods:

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
# # %%

# %%
