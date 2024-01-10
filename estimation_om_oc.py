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

# methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994',
#            'Malm_1994', 'Chow_1996', 'Andrews_2000',
#            'Malm_2000', 'Maenhaut_2002', 'DeBell_2006',
#            'Hand_2011', 'Simon_2011']
methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994',
           'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006',
           'Hand_2011', 'Simon_2011']

columnName="Event_M"
for method in methods:
    resultNormal = estimation_om_oc(matrix.where(
        events[columnName] == 'no'), method=method, ssa_as_Na=False)
    resultEvent = estimation_om_oc(matrix.where(
        events[columnName].isin(["S", "SP", "SN", "DS"])), method=method,
        ssa_as_Na=False)
    resultAll = estimation_om_oc(matrix, method=method,
        ssa_as_Na=False)
    #print(result.summary())

    # matrix.where(events[columnName].isin(["S", "SP", "SN","DS"]))[
    #     "PM2.5"].plot(style='o')

    warm_season_index = matrix.index.where(matrix.index.month >= 9).dropna()
    result = estimation_om_oc(matrix.loc[warm_season_index])
    print(f"{method}")
    print("No events")
    print(resultNormal.summary())#.as_latex())
    print(method, '&', f'{resultNormal.params[1]:.4g}', '&',
         f'{resultNormal.bse[1]:.2g}', '&',
         f'{resultNormal.pvalues[1]:.3g}', '\\\\')
    print("Events")
    print(resultEvent.summary())#.as_latex())
    print(method, '&', f'{resultEvent.params[1]:.4g}', '&',
         f'{resultEvent.bse[1]:.2g}', '&',
         f'{resultEvent.pvalues[1]:.3g}', '\\\\')
    # print("All")
    # print(resultAll.summary())#.as_latex())
    # print(method, '&', f'{resultAll.params[1]:.4g}', '&',
    #       f'{resultAll.bse[1]:.2g}', '&',
    #       f'{resultAll.pvalues[1]:.3g}', '\\\\')