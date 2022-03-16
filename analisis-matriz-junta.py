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
#Carga datos

df = pd.read_excel("Matriz.xlsx", sheet_name="valores", skiprows=[1,2], decimal=',').replace('#VALUE!', np.nan).replace("#DIV/0!", np.nan)
#df = df.where(df > 0)
df['Date'] = pd.to_datetime(df['Date'], format="%Y-%m-%d %H:%M:%S")


# +
#Carga meteorologia
obsbaires = pd.read_csv('../../doctorado/paper_laura/ARCAL_ISH/obsbaires.csv', delimiter=';') # obsbaires2 es abrir y volver a cerrar en excel para que ponga bien los delimiters, si no mezcla ; y ,
obsbaires = obsbaires[(obsbaires['date[yyyymmddHHMM]'] >= 201904031200)]
fechas = df['Date']
# Reemplazo los valores que son nan
obsbaires['wdir'] = obsbaires['wdir'].where(obsbaires['wdir'] < 999)
obsbaires['wspd[m/s]'] = obsbaires['wspd[m/s]'].where(obsbaires['wspd[m/s]'] < 999)
obsbaires['clht[km]'] = obsbaires['clht[km]'].where(obsbaires['clht[km]'] < 99)
obsbaires['dptp[C]'] = obsbaires['dptp[C]'].where(obsbaires['dptp[C]'] < 999)
obsbaires['slvp[hPa]'] = obsbaires['slvp[hPa]'].where(obsbaires['slvp[hPa]'] < 9999)
obsbaires['press[hPa]'] = obsbaires['press[hPa]'].where(obsbaires['press[hPa]'] < 9999)
obsbaires['prcp[mm]'] = obsbaires['prcp[mm]'].where(obsbaires['prcp[mm]'] < 999)
obsbaires['sky[octas]'] = obsbaires['sky[octas]'].where(obsbaires['sky[octas]'] > 99)
obsbaires['skyOpaque[octas]'] = obsbaires['skyOpaque[octas]'].where(obsbaires['skyOpaque[octas]'] < 99)
obsbaires['tmpd[C]'] = obsbaires['tmpd[C]'].where(obsbaires['tmpd[C]'] < 99)

# Calculo RH como rh=100*(EXP((17.625*d)/(243.04+d))/EXP((17.625*t)/(243.04+t)));
# donde d=dewpoint (dptp) y t=dry temp (tmpd), aprox de August-Roche-Magnus
d = obsbaires['dptp[C]']
t = obsbaires['tmpd[C]']
obsbaires['RH[%]'] = 100 * (np.exp((17.625 * d)/(243.04 + d)) / np.exp((17.625 * t)/(243.04 + t)))
# 

#Promedio fechas que están dobles
obsbaires = obsbaires.groupby('date[yyyymmddHHMM]').mean().reset_index()



# Parsing de fecha
datetime = pd.to_datetime(obsbaires['date[yyyymmddHHMM]'], format='%Y%m%d%H%M')
obsbaires['date'] = datetime
# Reordeno columnas-
cols = obsbaires.columns.tolist()
cols = cols[-1:] + cols[:-1] # Mando la última columna (date) al principio
obsbaires = obsbaires[cols]
#display(obsbaires)
obsbaires = obsbaires.set_index('date')

#display(obsbaires)

obsbaires = obsbaires.groupby(level=0).sum()
#display(obsbaires)

obsbaires24hs_mean = obsbaires.resample('24H', offset='12h').mean() # Resampleo para tener promedios diarios

#display(obsbaires24hs_mean)
obsbaires24hs_sum = obsbaires.resample('24H', offset='12h').sum() # Resampleo para tener promedios diarios

obsbaires_comb = pd.concat([obsbaires24hs_mean[['wdir', 'wspd[m/s]', 'clht[km]', 'hzvs[km]', 'tmpd[C]',
                                          'slvp[hPa]', 'press[hPa]', 'sky[octas]', 'skyOpaque[octas]', 'RH[%]']],
                         obsbaires24hs_sum[['prcp[mm]', 'prcpPeriod[hours]']]], axis=1)

#display(obsbaires_comb)
obsbaires_comb.index = obsbaires_comb.index.normalize()
obsbaires_comb = (obsbaires_comb.reindex(index = fechas.to_list()))
obsbaires_comb.to_excel('obsbaires_meteo2.xlsx')
# 
# #with pd.
# #display(fechas['obsbaires'])

# + tags=[] jupyter={"source_hidden": true}
matriz_de_corr = corr_matrix(df)
matriz_de_corr = matriz_de_corr.copy().round(decimals=2)

matriz_de_corr.to_csv('matriz_de_corr.csv', float_format='%g')

# + jupyter={"source_hidden": true} tags=[]
df['PM2.5'].plot(style='.-')
df['C Orgánico'].plot(style='.-')
df['C Total'].plot(style='.-')
plt.show()

plt.scatter(df['PM2.5'], df['C Orgánico'])
plt.xlabel('PM2.5')
plt.ylabel('OC')
print(df.keys())

# + tags=[]
keys = df.keys()

for i in range(1, len(keys)-1):
    for j in range(i+1, len(keys)-1):
        plt.scatter(df[keys[i]], df[keys[j]])
        plt.grid()
        plt.xlabel(keys[i])
        plt.ylabel(keys[j])
        plt.savefig(f'corr_nuevas/{keys[i]}_vs_{keys[j]}.png')
        plt.show()

# + jupyter={"outputs_hidden": true} tags=[]
for i in df.keys():
    df[i].plot(label=i)
    plt.legend()
    plt.show()

# +
df['Sr'] = 0
df['S'] = df['SO4']
methods = ['Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002',
           'DeBell_2006',
           'Hand_2011', 'Hand_2011_mod']

reconstruccion_masica = {}

for i in methods:
    reconstruccion_masica[i] = mass_closure(df, i)

for m in methods:
    escapan = 0
    escapan_defecto = 0
    escapan_exceso = 0
    cumplen = 0
    son_nan = 0
    criterio = 0.2
    for i in range(len(df)):
        value = ((mass_closure(data_df=df, equation=m)[0]) / df['PM2.5'])[i]
        if (1-criterio) < value < (1+criterio):
#           print('Filtro', i, value)
            cumplen += 1
        elif np.isnan(value):
#           print('Filtro', i, 'is nan')
            son_nan +=1
        else:
#            print('Filtro', i, value, 'escapa')
            escapan += 1
            if (1+criterio) < value:
                escapan_exceso += 1
            else:
                escapan_defecto += 1
    print(m + '\ncumplen =', cumplen, '\nescapan =', escapan,
          '\n\tpor defecto =', escapan_defecto, '\n\tpor exceso =', escapan_exceso,
          '\nson nan =', son_nan, '\n')

for i in methods:
    print('\n\n'+ str(i))
    mass = mass_closure(data_df=df, equation=i)
    print('Total explained mass = ', ((mass[0] / df['PM2.5']).mean()* 100).round(1), '%')
    for key in mass[1].keys():
        print(key, '=', ((mass[1][key] / df['PM2.5']).mean() * 100).round(1), '%') 

plt.figure(figsize=[20,10])
((reconstruccion_masica['Hand_2011'][0])/df['PM2.5']).plot(style='o-', label='Hand')
((reconstruccion_masica['Hand_2011_mod'][0])/df['PM2.5']).plot(style='o-', label='Hand_mod')
((reconstruccion_masica['DeBell_2006'][0])/df['PM2.5']).plot(style='o-', label='DeBelle')
plt.plot([0,120],[0.8,0.8], ':')
plt.plot([0,120],[1.2,1.2], ':')
plt.legend()
plt.show()
plt.plot(df['PM2.5'], reconstruccion_masica['Hand_2011'][0]/df['PM2.5'], 'o')
plt.xlabel('PM2.5')
plt.ylabel('Fracción explicada')
plt.show()

