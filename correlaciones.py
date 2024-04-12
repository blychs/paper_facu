import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from load_data import load_data
import seaborn as sn
from scipy import stats

%matplotlib
plt.style.use('seaborn-v0_8-paper')
matrix, unc, meteo, gases, events = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                              'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                              'BA_events_testM.xlsx')

print(matrix.keys())

matrix = matrix.drop(["ECPk1 C", "ECPk2 C", "ECPk3 C", "ECPk4 C", "ECPk5 C", "ECPk6 C"], axis=1)
matrix = matrix.drop(["OCPk1 C", "OCPk2 C", "OCPk3 C", "OCPk4 C", "OCPk5 C"], axis=1)
matrix = matrix.drop(["Pyrol C", "Na no sol", "Na total"], axis=1)
matrix = matrix.drop(["temp", "pres", "rh", "ws", "VentCoef"], axis=1)

print(matrix)

# filtered matrix, based on
# https://stackoverflow.com/questions/23199796/detect-and-exclude-outliers-in-a-pandas-dataframe

matrix = matrix[(np.abs(stats.zscore(df)) < 3).all(axis=1)]
keys = list(matrix.keys())
for i in range(len(keys) - 1):
    for j in range(i + 1, len(keys)):
        fig, ax = plt.subplots()
        ax.plot(matrix[keys[i]], matrix[keys[j]], 'o')
        fig.savefig(f"images/correlations_no_outliers/{keys[i]}_{keys[j]}_3sigma.png")
        plt.close()
    