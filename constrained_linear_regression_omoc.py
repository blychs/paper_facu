#%%
from constrained_linear_regression import ConstrainedLinearRegression
from sklearn.datasets import load_boston
from sklearn.linear_model import LinearRegression
# from sklearn.metrics import Regre
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
import statsmodels
import statsmodels.api as sm
import pandas as pd
import numpy as np
from sklearn import datasets, linear_model
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
from scipy import stats

plt.style.use('seaborn-v0_8-paper')
matrix, unc, meteo, gases, events = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                              'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                              'BA_events_testM.xlsx')
datesdrop=['2019-05-24','2019-05-27','2019-05-30','2019-06-02', '2020-03-01','2020-01-31']
matrix=matrix.drop(datesdrop,axis=0)
events=events.drop(datesdrop,axis=0)

methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994',
           'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006',
           'Hand_2011','Simon_2011']

event_columnname="Event_M"
event_labels= ["S", "SP", "SN","SL","SC"]
omoc_noevent=[]
omoc_event=[]
omoc_all=[]
method='Simon_2011'
ssa_as_Na=False
display_latex=False
concentration_matrix = matrix.copy()
if "C Elemental" in concentration_matrix:
    concentration_matrix["EC"] = concentration_matrix["C Elemental"]
if "C Orgánico" in concentration_matrix:
    concentration_matrix["OC"] = concentration_matrix['C Orgánico']
if "Si" not in concentration_matrix:
    concentration_matrix['Si'] = 2.4729 * concentration_matrix['Al']
if "(NH4)2SO4" not in concentration_matrix:
    # Assume all SO4 is (NH4)2SO4
    concentration_matrix["(NH4)2SO4"] = (132.14 / 96.06 ) * concentration_matrix["SO4"]
if "NH4NO3" not in concentration_matrix:
    # Assume all NO3 is NH4NO3
    concentration_matrix["NH4NO3"] =  (80.043 / 62.004 ) * concentration_matrix["NO3"]
    

concentration_matrix = concentration_matrix.dropna(subset=["(NH4)2SO4", "NH4NO3", "Si", "Ca", 
                                                    "Fe", "Ti", "EC", "Cl","K", "Fe", 
                                                    "PM2.5", "OC"], axis=0)

soil = (3.48 * concentration_matrix["Si"] + 1.63 * concentration_matrix["Ca"] +
        2.42 * concentration_matrix["Fe"] + 1.94 * concentration_matrix["Ti"])

if not ssa_as_Na:
    salt = 1.8 * concentration_matrix["Cl"]
else:
    salt = 2.54 * concentration_matrix["Na sol"]
intercept_base = (concentration_matrix["EC"] + salt +
                    1.2 * (concentration_matrix["K"] - 0.6 * concentration_matrix["Fe"]))

#    print(concentration_matrix)
model = ConstrainedLinearRegression()
min_coef = np.array([1,0.59,-0.9,0.41])
# min_coef = np.array([-np.inf,-np.inf,-np.inf,-np.inf])
#min_coef = np.array([-np.inf,1,1,1])

max_coef = np.array([3.8,1.53,1.35,1.63])
# max_coef = np.array([np.inf,np.inf,np.inf,np.inf])
#max_coef = np.array([np.inf,1,1,1])


print("Events")        
y = (concentration_matrix['PM2.5'] - intercept_base).where(events[event_columnname].isin(event_labels)).dropna().values
X = np.column_stack((concentration_matrix["OC"].where(events[event_columnname].isin(event_labels)).dropna().values,
                    concentration_matrix["(NH4)2SO4"].where(events[event_columnname].isin(event_labels)).dropna().values,
                    concentration_matrix["NH4NO3"].where(events[event_columnname].isin(event_labels)).dropna().values,
                    soil.where(events[event_columnname].isin(event_labels)).dropna().values))
model.fit(X, y, max_coef=max_coef, min_coef=min_coef)
# print(model.intercept_)
print(model.coef_)
betas_event=model.coef_
print_stats(model, X, y)

print("No Events")    
y = (concentration_matrix['PM2.5'] - intercept_base).where(events[event_columnname] == 'no').dropna().values
X = np.column_stack((concentration_matrix["OC"].where(events[event_columnname] == 'no').dropna().values,
                    concentration_matrix["(NH4)2SO4"].where(events[event_columnname] == 'no').dropna().values,
                    concentration_matrix["NH4NO3"].where(events[event_columnname] == 'no').dropna().values,
                    soil.where(events[event_columnname] == 'no').dropna().values))
model.fit(X, y, max_coef=max_coef, min_coef=min_coef)
# print(model.intercept_)
print(model.coef_)
betas_noevent=model.coef_
print_stats(model, X, y)

print("All together")    
y = (concentration_matrix['PM2.5'] - intercept_base).values
X = np.column_stack((concentration_matrix["OC"].values,
                    concentration_matrix["(NH4)2SO4"].values,
                    concentration_matrix["NH4NO3"].values,
                    soil.values))
model.fit(X, y, max_coef=max_coef, min_coef=min_coef)
print(model.coef_)
betas_all=model.coef_
print_stats(model, X, y)


# %%

