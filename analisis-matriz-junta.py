# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
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
#Carga datos

df_conc = pd.read_excel("PMF_BA_CONCyUNC.xlsx", decimal=',').replace(-999., np.nan)
#df = df.where(df > 0)
df_conc['date'] = pd.to_datetime(df_conc['date'], format="%Y-%m-%d %H:%M:%S")

#df_unc = pd.read_excel("PMF_BA_completo.xlsx", sheet_name="UNC", decimal=',').replace(-999., np.nan)
#df_unc['date'] = pd.to_datetime(df_unc['date'], format="%Y-%m-%d %H:%M:%S")

df_conc = df_conc.rename(columns={'PM2,5': 'PM2.5'})
#df_unc = df_unc.rename(columns={'PM2,5': 'PM2.5'})

#(df_conc['uC Orgánico']/df_conc['C Orgánico'] * 100).plot(label='eOC')
#(df_conc['uC Elemental']/df_conc['C Elemental'] * 100).plot(label='eEC')
#(df_conc['uNa']/df_conc['Na'] * 100).plot(label='eNa')

#plt.legend()
#plt.show()

#exclude = ['date', 'C Total', 'uC Elemental', 'uC Orgnico', 'Na total']
#df_conc.loc[:, df_conc.columns.difference(exclude)].plot()

#display(df_unc)

#df_normal = df.where(df['E. Regional binario'] == 'no')
#df_evento = df.where(df['E. Regional binario'] == 'si')
display(df_conc)

print(df_conc.keys().to_list())


# + tags=[]
values_for_corr = ['PM2.5', 'Cl', 'NO3', 'SO4', 'Na', 'NH4', 'C Orgánico',
                   'C Elemental', 'C Total', 'Na total', 'Mg', 'Al', 'K',
                   'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                   'Zn', 'As', 'Se', 'Mo', 'Ag', 'Cd', 'Sb', 'Ba', 'Pb']


matriz_de_corr_p = corr_matrix(df_conc[values_for_corr], method='pearson')
matriz_de_corr_p = matriz_de_corr_p.copy().round(decimals=2)
matriz_de_corr_p.to_csv('matriz_de_corr_pearson.csv', float_format='%g')

matriz_de_corr_s = corr_matrix(df_conc[values_for_corr], method='spearman')
matriz_de_corr_s = matriz_de_corr_s.copy().round(decimals=2)
matriz_de_corr_s.to_csv('matriz_de_corr_spearman.csv', float_format='%g')

matriz_de_corr_k = corr_matrix(df_conc[values_for_corr], method='kendall')
matriz_de_corr_k = matriz_de_corr_k.copy().round(decimals=2)
matriz_de_corr_k.to_csv('matriz_de_corr_kendall.csv', float_format='%g')


# + tags=[]
#Grafico de corrs
# keys = df.keys()
# 
# for i in range(4, len(keys)-1):
#     for j in range(i+1, len(keys)-1):
#         plt.scatter(df_normal[keys[i]], df_normal[keys[j]])
#         plt.scatter(df_evento[keys[i]], df_evento[keys[j]])
#         plt.grid()
#         plt.xlabel(keys[i])
#         plt.ylabel(keys[j])
#         plt.savefig(f'corr_nuevas/{keys[i]}_vs_{keys[j]}.png')
#         plt.show()

# +
def porcentaje(f, PM):
    return f/PM * 100

def err_porcentaje(f, PM, err_f, err_PM):
    return np.linalg.norm( [ 100/PM *err_f, f/(PM**2) * 100 * err_PM ], axis=0)


methods = ['Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002',
           'DeBell_2006',
           'Hand_2011', 'Hand_2011_mod']

for m in methods:
    closure = mass_closure(df_conc, equation=m)
    f = closure[0]
    PM = df_conc['PM2.5']
    err_f = closure[2]
    err_PM = df_conc['uPM2,5']
    
    fig, ax1 = plt.subplots(figsize=(15,7.5))
    color = 'tab:blue'
    ax1.set_xlabel('Date', size=14)
    ax1.set_ylabel('Reconstructed mass (%)', color=color, size=14)
    ax1.tick_params(axis='y', labelcolor=color, labelsize=12)
    ax1.tick_params(axis='x', labelsize=12)
    ax1.errorbar(df_conc['date'], porcentaje(f, PM), err_porcentaje(f, PM, err_f, err_PM),
                 capsize=3, ecolor='r', fmt='o-', label='Percentage reconstructed')
    ax1.set_ylim(ymin=0)
    #ax1.axhline(y=80, color='k', linestyle=':')
    #ax1.axhline(y=120, color='k', linestyle=':')
    ax1.axhspan(80, 120, color='y', alpha=0.2)
    ax1.legend(loc=2, frameon=False, fontsize=12)
    ax1.set_title(f'Method = {m}')
    # Adding values
    
    ax2 = ax1.twinx()
    color = 'tab:red'
    ax2.set_ylabel('Absolute concentration (μg/m$^3$)', color=color, size=14)
    ax2.tick_params(axis='y', labelcolor=color, labelsize=12)
    ax2.errorbar(df_conc['date'], df_conc['PM2.5'], df_conc['uPM2,5'],
                 capsize=4, ecolor='r', color='darkred', fmt='o-', label='PM2.5 measured')
    ax2.errorbar(df_conc['date'], closure[0], closure[2],
                 capsize=4, ecolor='red', color='fuchsia', fmt='o:',
                 label='Reconstructed mass')
    ax2.set_ylim(ymin=0, ymax=180)
    ax2.legend(frameon=False, fontsize=12)
    
    plt.savefig(f'{m}.png')
    
    plt.show()


#plt.figure(figsize=(20, 10))
#plt.errorbar(df_conc['date'], porcentaje(f, PM), err_porcentaje(f, PM, err_f, err_PM), capsize=3, ecolor='r', fmt='o-')
#plt.plot(df_conc['date'], df_conc['PM2.5'])
#plt.plot(df_conc['date'], closure[0], ':')
#plt.axhline(y=80, color='k', linestyle=':')
#plt.axhline(y=120, color='k', linestyle=':')
#plt.ylim(ymin=0)

# + tags=[]
dfresults=pd.DataFrame(columns=["method", "explained", "too_much", "too_little", "nan"])

for m in methods:
    dresults = {"method": m, "explained": 0, "too_much": 0, "too_little":0, "nan": 0}
    criterio = 20
    cierre = mass_closure(data_df=df_conc, equation=m)
    cierre_por = porcentaje(cierre[0], df_conc['PM2.5'])
    error_por = err_porcentaje(cierre[0], df_conc['PM2.5'], cierre[2], df_conc['uPM2,5'])
    
    for i in range(len(cierre_por)):
        value = cierre_por[i]
        uvalue = error_por[i]
        if (100 - criterio) < value < (100 + criterio):
            dresults["explained"] += 1
        elif np.isnan(value):
            dresults["nan"] +=1
        elif value < 100 - criterio:
            value += uvalue
            if value < 100 - criterio:
                dresults["too_little"] +=1
            else:
                dresults["explained"] += 1
        elif value > (100 + criterio):
            value -= uvalue
            if value > 100 + criterio:
                dresults["too_much"] += 1
            else:
                dresults["explained"] += 1
    
    dfresults = pd.concat([dfresults, pd.DataFrame([dresults])], ignore_index=True)
dfresults.set_index("method", inplace=True)
display(dfresults)


# +
# Percentage of measured mass by category

dcolors = {'inorganic_ions': 'b', 'organic_mass': 'r', 'elemental_C': 'c', 'trace_elements': 'm', 
           'geological_minerals': 'tab:brown', 'salt': 'silver', 'others':'violet'}

for m in methods:
    closure = mass_closure(df_conc, equation=m)
    
    PM = df_conc['PM2.5']
    err_PM = df_conc['uPM2,5']
    fig, ax1 = plt.subplots(figsize=(15,7.5))
    fig, ax2 = plt.subplots(figsize=(15,7.5))
    ax2.errorbar(df_conc['date'],df_conc['PM2.5'],
                  df_conc['uPM2,5'], color='k', ecolor='k',
                  capsize=3, fmt='o-', label='PM2.5')
    for j in list(closure[1].keys())[:-1]:
        f = closure[1][j]
        err_f = closure[3][f'u{j}']
        color = 'tab:blue'
        ax1.set_xlabel('Date', size=14)
        ax1.set_ylabel('Mass (%)', color=color, size=14)
        ax1.tick_params(axis='y', labelcolor=color, labelsize=12)
        ax1.tick_params(axis='x', labelsize=12)
        ax1.errorbar(df_conc['date'], porcentaje(f, PM),
                     err_porcentaje(f, PM, err_f, err_PM), color=dcolors[j], ecolor=dcolors[j],
                     capsize=3, fmt='o-', label=j)
        ax1.set_ylim(ymin=0, ymax=100)
        #ax1.axhline(y=80, color='k', linestyle=':')
        #ax1.axhline(y=120, color='k', linestyle=':')
        ax1.legend(loc=2, frameon=False, fontsize=12)
        ax1.set_title(f'Method = {m}')
        # Adding values
        ax2.set_xlabel('Date', size=14)
        ax2.set_ylabel('Mass', color=color, size=14)
        ax2.tick_params(axis='y', labelcolor=color, labelsize=12)
        ax2.tick_params(axis='x', labelsize=12)
        ax2.errorbar(df_conc['date'],f,
                     err_f, color=dcolors[j], ecolor=dcolors[j],
                     capsize=3, fmt='o-', label=j)
        ax2.set_ylim(ymin=0)
        #ax1.axhline(y=80, color='k', linestyle=':')
        #ax1.axhline(y=120, color='k', linestyle=':')
        ax2.legend(loc=2, frameon=False, fontsize=12)
        ax2.set_title(f'Method = {m}, abs')
    
        
    
    plt.show()
# -

# Percentage of reconstructed mass
for m in methods:
    closure = mass_closure(df_conc, equation=m)
    
    PM = closure[0]
    err_PM = closure[2]
    fig, ax1 = plt.subplots(figsize=(15,7.5))
    fig, ax2 = plt.subplots(figsize=(15,7.5))
    for j in list(closure[1].keys())[:-1]:
        f = closure[1][j]
        err_f = closure[3][f'u{j}']
        color = 'tab:blue'
        ax1.set_xlabel('Date', size=14)
        ax1.set_ylabel('Reconstructed mass (%)', color=color, size=14)
        ax1.tick_params(axis='y', labelcolor=color, labelsize=12)
        ax1.tick_params(axis='x', labelsize=12)
        ax1.errorbar(df_conc['date'], porcentaje(f, PM),
                     err_porcentaje(f, PM, err_f, err_PM),
                     capsize=3, fmt='o-', label=j)
        ax1.set_ylim(ymin=0, ymax=100)
        #ax1.axhline(y=80, color='k', linestyle=':')
        #ax1.axhline(y=120, color='k', linestyle=':')
        ax1.legend(loc=2, frameon=False, fontsize=12)
        ax1.set_title(f'Method = {m}')
        # Adding values
        ax2.set_xlabel('Date', size=14)
        ax2.set_ylabel('Reconstructed mass', color=color, size=14)
        ax2.tick_params(axis='y', labelcolor=color, labelsize=12)
        ax2.tick_params(axis='x', labelsize=12)
        ax2.errorbar(df_conc['date'],f,
                     err_f,
                     capsize=3, fmt='o-', label=j)
        ax2.set_ylim(ymin=0)
        #ax1.axhline(y=80, color='k', linestyle=':')
        #ax1.axhline(y=120, color='k', linestyle=':')
        ax2.legend(loc=2, frameon=False, fontsize=12)
        ax2.set_title(f'Method = {m}, abs')
    
        
    
    plt.show()
