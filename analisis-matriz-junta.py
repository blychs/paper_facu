# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import numpy as np # librería de cálculo numérico
from scipy import stats # librería estadística
import matplotlib.pyplot as plt # para gráficos
import pandas as pd # para tablas (tipo spreadsheets)
from pandas.api.types import is_numeric_dtype # chequeos de tipo numérico
import xarray as xr # es como una extensión de Pandas para usar varias dimensiones
                    # y hacer cálculos con vectores

# Funciones que creé yo
# %run ./funciones.ipynb

# + tags=[]
df = pd.read_excel("Matriz.xlsx", sheet_name="valores", skiprows=[1,2], decimal=',').replace('#VALUE!', np.nan).replace("#DIV/0!", np.nan)
df = df.where(df > 0)


# + tags=[]
matriz_de_corr = corr_matrix(df)
matriz_de_corr = matriz_de_corr.copy().round(decimals=2)

matriz_de_corr.to_csv('matriz_de_corr.csv', float_format='%g')

# +
df['PM2.5'].plot(style='.-')
df['C Orgánico'].plot(style='.-')
df['C Total'].plot(style='.-')
plt.show()

plt.scatter(df['PM2.5'], df['C Orgánico'])
plt.xlabel('PM2.5')
plt.ylabel('OC')

# + jupyter={"outputs_hidden": true} tags=[]
keys = df.keys()

for i in range(1, len(keys)-1):
    for j in range(i+1, len(keys)-1):
        plt.scatter(df[keys[i]], df[keys[j]])
        plt.grid()
        plt.xlabel(keys[i])
        plt.ylabel(keys[j])
        plt.savefig(f'corr_nuevas/{keys[i]}_vs_{keys[j]}.png')
        plt.show()

# +
methods = [#'Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           #'Malm_2000', #'Maenhaut_2002',
           #'DeBell_2006',
           'Hand_2011', 'Hand_2011_mod']

reconstruccion_masica = {}

for i in methods:
    reconstruccion_masica[i] = mass_closure(df, i)

plt.figure(figsize=[20,10])
((reconstruccion_masica['Hand_2011'][0])/df['PM2.5']).plot(style='o-', label='Hand')
((reconstruccion_masica['Hand_2011_mod'][0])/df['PM2.5']).plot(style='o-', label='Hand_mod')
plt.plot([0,120],[0.8,0.8], ':')
plt.plot([0,120],[1.2,1.2], ':')
plt.legend()



