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
# import xarray as xr # es como una extensión de Pandas para usar varias dimensiones
                    # y hacer cálculos con vectores

# Funciones que creé yo
# %run ./funciones.ipynb

# +
df_raw = pd.read_excel('20210511_matriz_36_fin.xlsx', index_col=0, skiprows=[1,2]) # Matriz CSV curada

############# Remove outliers ############
# Elimino outliers de la matriz definido como aquello que está a más de 3 sigmas 
# (falta pasar a log)
df_log = np.log(df_raw)
#df = df_raw.mask(df_raw.sub(df_raw.mean()).div(df_raw.std()).abs().gt(3)) 
df = df_raw.mask(df_log.sub(df_log.mean()).div(df_log.std()).abs().gt(3)) 

df_tecnicas = pd.read_csv('ArcalMetalesAnalisis.csv', index_col=0, skiprows=[1,2])


# + tags=[]
print(df['C Orgánico'])

# + jupyter={"source_hidden": true} tags=[]
## Ploteo de elementos claves que queramos ver para comparar técnicas
## (modificar lista_elementos si se quiere cambiar cuáles)
###########
#lista_elementos = ['Ti', 'Mn', 'Ni', 'Zn', 'Pb']
#for i in lista_elementos:
#    ax = comparacion_tecnicas(df_tecnicas, i)
#    plt.show()

# + tags=[]
# CÁLCULO DE MATRIZ DE CORRELACIONES
df.corr().to_csv('correlacion_todo.csv')
df.plot(label='')

# + tags=[] jupyter={"outputs_hidden": true}
# CÁLCULO DE LOGNORMALIDAD Y GRÁFICOS DE HISTOGRAMAS CLAVE

#Histogramas de la versión lineal
df.plot(y=['C Orgánico', 'C Elemental', 'C Total', 'PM2.5'], ylim=(0,30))
df.hist('C Orgánico', bins=20)
df.hist('C Elemental', bins=20)
df.hist('C Total', bins=20)
plt.show()

# Histogramas logarítmicos seleccionados
lista_de_elementos = ['PM2.5' ,'SO4', 'Pb', 'Ca']
for i in lista_de_elementos:
    log_i = np.log(df[i])
    log_i.hist(bins=10)
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
            if p_is_lognormal < 0.05:
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
# -

df.plot.scatter('PM2.5', 'C Total')
plt.show()
df.plot.scatter('PM2.5', 'C Elemental')
plt.show()
df.plot.scatter('PM2.5', 'C Orgánico')
plt.show()
df.plot.scatter('PM2.5', 'Na')
plt.show()
df.plot.scatter('PM2.5', 'Cl')
plt.show()

# +
#Ploteo de carbonosas y SO4

(df['C Elemental'] * 8).plot(label='Elemental * 8')
(df['PM2.5']).plot(label='PM2.5')
plt.legend()
plt.show()
(df['C Orgánico']).plot(label='Orgánico')
(df['PM2.5']).plot(label='PM2.5')
plt.legend()
plt.show()
(df['SO4'] * 8).plot(label='SO4 * 8')
(df['C Orgánico']).plot(label='C Org')
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

#for k in methods:
#    print('\n', k)
#    mass_closure_fraction = (mass_closure(data_df=df, equation=k)[0] / df['PM2.5'])
#    print(mass_closure_fraction)
    #for i in range(0, len(mass_closure_fraction)):
    #    print(mass_closure_fraction[0])
         #print('Filtro ' + str(i + 1) + ' = ' + str(mass_closure_fraction[i]))


        
for i in methods:
    print('\n\n',i)
    print(str((mass_closure(data_df=df, equation=i)[0] / df['PM2.5']).mean()))
    print('inorganic_ions ', (mass_closure(data_df=df, equation=i)[1] / df['PM2.5']).mean())
    print('organic_mass ', (mass_closure(data_df=df, equation=i)[2] / df['PM2.5']).mean())
    print('EC ', (mass_closure(data_df=df, equation=i)[3] / df['PM2.5']).mean())
    print('geological_minerals ', (mass_closure(data_df=df, equation=i)[4] / df['PM2.5']).mean())
    print('seasalt ', (mass_closure(data_df=df, equation=i)[5] / df['PM2.5']).mean())
    print('trace_elements ', (mass_closure(data_df=df, equation=i)[6] / df['PM2.5']).mean())
    print('others ', (mass_closure(data_df=df, equation=i)[7] / df['PM2.5']).mean())
    print('unexplained ', (mass_closure(data_df=df, equation=i)[8] / df['PM2.5']).mean())
    

# closure, inorganic_ions, organic_mass, elemental_C,
#geological_minerals, salt, trace_elements, others, unexplained =




# +
escapan = 0
cumplen = 0
son_nan = 0
criterio = 0.2
for i in range(38,120):
    value = ((mass_closure(data_df=df, equation='Hand_2011')[0]) / df['PM2.5'])[i]
    if (1-criterio) < value < (1+criterio):
#        print('Filtro', i, value)
        cumplen += 1
    elif np.isnan(value):
#        print('Filtro', i, 'is nan')
        son_nan +=1
    else:
#        print('Filtro', i, value, 'escapa')
        escapan += 1

print('Hand 2011\ncumplen =', cumplen, '\nescapan =', escapan, '\nson nan =', son_nan)

escapan = 0
cumplen = 0
son_nan = 0
criterio = 0.2
for i in range(38,120):
    value = ((mass_closure(data_df=df, equation='Hand_2011_mod')[0]) / df['PM2.5'])[i]
    if (1-criterio) < value < (1+criterio):
#        print('Filtro', i, value)
        cumplen += 1
    elif np.isnan(value):
#        print('Filtro', i, 'is nan')
        son_nan +=1
    else:
#        print('Filtro', i, value, 'escapa')
        escapan += 1

print('\n\nHand 2011 mod\ncumplen =', cumplen, '\nescapan =', escapan, '\nson nan =', son_nan)


style='o-'
ax, figure = plt.subplots(figsize=(20,10))
(mass_closure(data_df=df, equation='Hand_2011')[0] / df['PM2.5']).plot(style=style, label='Hand')
((mass_closure(data_df=df, equation='Hand_2011_mod')[0] / df['PM2.5']).plot(style=style, color='r', label='Hand modif'))
#((mass_closure(data_df=df, equation='DeBell_2006')[0] / df['PM2.5']).plot(style='o-', color='g', label='DeBell'))
#((mass_closure(data_df=df, equation='Malm_2000')[0] / df['PM2.5']).plot(style='o-', color='k', label='Malm'))
plt.legend()

plt.plot([38,120], [1+criterio, 1+criterio])
plt.plot([38,120], [1-criterio, 1-criterio])
plt.show()

# -

value = ((mass_closure(data_df=df, equation='Hand_2011')[0]) / df['PM2.5'])[41]
print(np.isnan(value))

# + tags=[] jupyter={"outputs_hidden": true}
# Tortas
methods = ['Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011']

labels = 'II', 'OM', 'EC','GM', 'Sea S', 'TE', 'O', 'NC'

#sizes = (mass_closure(data_df=df.iloc[1], equation=k)[1:])
for k in methods:
    print(k)
    for j in range(len(df)):
        print(j)
        sizes = (mass_closure(data_df=df.iloc[j], equation=k)[1:])
   # print(len(sizes))
    #explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')
        try:
            fig1, ax1 = plt.subplots()
            ax1.pie(sizes, labels=labels)
            ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
            plt.savefig('imagenes/tortas/torta_filtro_' + 
                        str(df.index[j]) + '_metodo_' + k + '.eps')
            plt.show()
        except:
            print('no se pudo hacer', j)

# + tags=[]
# Correlaciones de todo con todo
graph_all_corr(df, 'correlaciones/')
# +
df.keys()
lista_keys = ['Cl', 'NO3', 'SO4', 'Na', 'NH4', 'C Orgánico', 'C Elemental',
              'C Total', 'S', 'Ca', 'V', 'Cr', 'Fe', 'Cu', 'As', 'Se', 'Sr',
              'Pb', 'Zn', 'Mn', 'Mo', 'Ni', 'Ti', 'Sb', 'PM2.5']

for i in range(len(lista_keys)):
    for j in range(i+1, len(lista_keys)):
        df_i_norm = df[lista_keys[i]] / df[lista_keys[i]].mean()
        df_j_norm = df[lista_keys[j]] / df[lista_keys[j]].mean()
        df_i_norm.plot(label=lista_keys[i])
        df_j_norm.plot(label=lista_keys[j])
        plt.legend()
        plt.show()


# +
df.keys()
lista_keys = ['Cl', 'NO3', 'SO4', 'Na', 'NH4', 'C Orgánico', 'C Elemental',
              'C Total', 'S', 'Ca', 'V', 'Cr', 'Fe', 'Cu', 'As', 'Se', 'Sr',
              'Pb', 'Zn', 'Mn', 'Mo', 'Ni', 'Ti', 'Sb', 'PM2.5']

for i in range(len(lista_keys)):
    df_i_norm = df[lista_keys[i]] / df[lista_keys[i]].mean()
    df_PM25_norm = df['PM2.5'] / df['PM2.5'].mean()
    df_i_norm.plot(label=lista_keys[i])
    df_PM25_norm.plot(label='PM2.5')
    plt.legend()
    plt.show()
