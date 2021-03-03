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
#     display_name: Python 3
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

# +
df_raw = pd.read_csv('20210218_matriz_frankestein.csv', index_col=0, skiprows=[1,2])[33:62] # Matriz CSV curada

############# Remove outliers ############
# Elimino outliers de la matriz definido como aquello que está a más de 3 sigmas 
# (falta pasar a log)
df = df_raw.mask(df_raw.sub(df_raw.mean()).div(df_raw.std()).abs().gt(3)) 

df_tecnicas = pd.read_csv('ArcalMetalesAnalisis.csv', index_col=0, skiprows=[1,2])


# +
# Ploteo de elementos claves que queramos ver para comparar técnicas
# (modificar lista_elementos si se quiere cambiar cuáles)
###########
# lista_elementos = ['Ti', 'Mn', 'Ni', 'Zn', 'Pb']
# for i in lista_elementos:
#     ax = comparacion_tecnicas(df_tecnicas, i)
#   plt.show()
# -

# CÁLCULO DE MATRIZ DE CORRELACIONES
df.corr().to_csv('correlacion_todo.csv')

# +
# CÁLCULO DE LOGNORMALIDAD Y GRÁFICOS DE HISTOGRAMAS CLAVE

#Histogramas de la versión lineal
df.plot(y=['C Orgánico', 'C Elemental', 'C Total', 'PM2.5'], ylim=(0,20))
df.hist('C Orgánico', bins=20)
df.hist('C Elemental', bins=20)
df.hist('C Total', bins=20)
plt.show()

# Histogramas logarítmicos seleccionados
lista_de_elementos = ['PM2.5' ,'SO4', 'Pb', 'Ca']
for i in lista_de_elementos:
    log_i = np.log(df[i])
    log_i.hist(bins=20)
    plt.title('ln(' + i + ')')
    plt.show()


# +
not_lognormal = []
lognormality_nan = []

### This function tests the null hypothesis that a sample comes from a normal distribution.
### It is based onD’Agostino and Pearson’s [1], [2] test that combines skew and kurtosis 
### to produce an omnibus test of normality.

with open('lognormalidad.txt', 'w') as lognorm:
    for i in df:
        if not df[i].isna().all() and i != 'Se.1': 
            p_is_lognormal = test_lognormality(df[i])
            print(i, str(p_is_lognormal))
            if p_is_lognormal < 0.01:
                not_lognormal.append(i)
            elif p_is_lognormal != p_is_lognormal: #Check if it's np.nan
                lognormality_nan.append(i)
        lognorm.write(i +' '+ str(p_is_lognormal) + '\n')
    lognorm.write('No son lognormales ' + str(not_lognormal) + '\n')
    lognorm.write('No hay valores suficientes para evaluar' + str(lognormality_nan) + '\n')
print('\n\n')
print('No son lognormales ' + str(not_lognormal))

# +
#for i in df:
#    if is_numeric_dtype(df[i]):
#        #plt.figure(figsize=(20,10))
#        df.plot(y=i, label=i, xlim=(0, 96))
#        plt.show()

# +
#Ploteo de carbonosas y SO4

(df['C Elemental'] * 8).plot(label='Elemental * 16')
(df['PM2.5']).plot(label='PM2.5', xlim=(0,70))
plt.legend()
plt.show()
(df['C Orgánico']).plot(label='Orgánico')
(df['PM2.5']).plot(label='PM2.5', xlim=(0,70))
plt.legend()
plt.show()
(df['SO4'] * 8).plot(label='SO4 * 8')
(df['C Orgánico']).plot(label='C Org', xlim=(0,70))
plt.legend()
plt.show()

# +
# Fracción de orgánico y de elemental
fraccion_org = (df['C Orgánico']*1.4) / df['PM2.5']
fraccion_elemental = df['C Elemental'] / df['PM2.5']
(fraccion_org).plot(label='C Orgánico/PM2.5')
plt.legend()
plt.show()
print(fraccion_org.mean())

#print(fraccion_elemental.mean())

# +

fraccion_NO3 = df['NO3'] / df['PM2.5']
(fraccion_NO3).plot(label='NO3/PM2.5')
plt.legend()
plt.show()
print(fraccion_NO3.mean())

# +
###### NO SÉ QUÉ ES ESTO, PASA POR NO COMENTAR#####
# for i in df:
#     if 'KED' not in i:
#         print(i)
#         print(df[i].mean()/df['PM2.5'].mean())
# 
# masa_explicada = 0
# for i in df:
#     if '.1' not in i:
#         masa_explicada += df[i].mean()
# masa_explicada = (masa_explicada - df['PM2.5'].mean() - df['C Total'].mean())/(df['PM2.5'].mean())
# print('masa explicada = ', str(masa_explicada))
##################################

# +
(df['C Total'] *1.4 * 2.5).plot(label='Orgánico')
(df['Na'] * 10).plot(label='Na')
(df['S']).plot(label='S', xlim=(0,40))
(df['NO3']).plot(label='NO$_3$')
(df['NH4']).plot(label='NH$_4$')
plt.legend()
plt.show()

df['S'].plot(label='S', xlim=(0,40))
df['SO4'].plot(label='SO$_4$')
plt.legend()
plt.show()

df['Ni.1'].plot(label='Ni.1')
df['Ni'].plot(label='Ni')
plt.legend()
plt.show()

df['Pb'].plot(label='Pb')
df['Sb'].plot(label='Sb')
plt.legend()
plt.show()
# -



# +

correlacion = corr_matrix(df)
correlacion.to_csv('correlacion.csv')
keys = correlacion.keys()

df[keys[0]]
# -

plt.plot((df['S'] - df['SO4'] * 0.33), df['V'], 'yo')
plt.show()
(df['V'] * 100).plot()
(df['S'] - df['SO4'] * 0.33).plot()

# +
methods = ['Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011']

for k in methods:
    print('\n', k)
    mass_closure_fraction = (mass_closure(data_df=df, equation=k)[0] / df['PM2.5'])
    print(mass_closure_fraction)
    #for i in range(0, len(mass_closure_fraction)):
    #    print(mass_closure_fraction[0])
         #print('Filtro ' + str(i + 1) + ' = ' + str(mass_closure_fraction[i]))


for i in methods:
    print('\n\n',i)
    print(str((mass_closure(data_df=df, equation=i)[0] / df['PM2.5']).mean()))
        


# +
# Tortas
methods = ['Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011']

labels = 'II', 'OM', 'EC','GM', 'Sea S', 'TE', 'O'

#sizes = (mass_closure(data_df=df.iloc[1], equation=k)[1:])
for k in methods:
    print(k)
    for j in range(len(df)):
        print(j)
        sizes = (mass_closure(data_df=df.iloc[j], equation=k)[1:])
   # print(len(sizes))
    #explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')

        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, labels=labels)
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        plt.savefig('')
        plt.show()
# -

# Correlaciones de todo con todo
graph_all_corr(df, 'correlaciones/')


