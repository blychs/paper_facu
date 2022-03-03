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
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd # para tablas (tipo spreadsheets)
from pandas.api.types import is_numeric_dtype # chequeos de tipo numérico
# import xarray as xr # es como una extensión de Pandas para usar varias dimensiones
                    # y hacer cálculos con vectores

# Funciones que creé yo
# %run ./funciones.ipynb
# #%matplotlib widget


df_raw = pd.read_excel('20210621_matriz_fechas.xlsx', index_col=0, skiprows=list(np.arange(1, 39)), usecols= "A,C:AL" ) # Matriz CSV curada

############# Remove outliers ############
# Elimino outliers de la matriz definido como aquello que está a más de 3 sigmas 
# (falta pasar a log)
df_log = np.log(df_raw)
#df = df_raw.mask(df_raw.sub(df_raw.mean()).div(df_raw.std()).abs().gt(3)) 
df = df_raw.mask(df_log.sub(df_log.mean()).div(df_log.std()).abs().gt(9)) 

df_tecnicas = pd.read_csv('ArcalMetalesAnalisis.csv', index_col=0, skiprows=[1,2])

# + tags=[]
display(df)

# + tags=[]
# CÁLCULO DE MATRIZ DE CORRELACIONES
df.corr().to_csv('correlacion_todo.csv')
df.plot(legend=True)
plt.show()

# + tags=[]
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


# + tags=[]
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

# + tags=[]
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

# + tags=[]
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

# + tags=[]
# Fracción de orgánico y de elemental
fraccion_org = (df['C Orgánico']*1.4) / df['PM2.5']
fraccion_elemental = df['C Elemental'] / df['PM2.5']
(fraccion_org).plot(label='C Orgánico/PM2.5')
plt.legend()
plt.show()
print(fraccion_org.mean())

#print(fraccion_elemental.mean())

# + tags=[]

fraccion_NO3 = df['NO3'] / df['PM2.5']
(fraccion_NO3).plot(label='NO3/PM2.5')
plt.legend()
plt.show()
print(fraccion_NO3.mean())

# + tags=[]
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

correlacion = corr_matrix(df)
correlacion.to_csv('correlacion.csv')
keys = correlacion.keys()

df[keys[0]]
# -

plt.plot((df['S'] - df['SO4'] * 0.33), df['V'], 'yo')
plt.show()
(df['V'] * 100).plot()
(df['S'] - df['SO4'] * 0.33).plot()

# + tags=[]
methods = ['Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011', 'Hand_2011_mod']

#for k in methods:
#    print('\n', k)
#    mass_closure_fraction = (mass_closure(data_df=df, equation=k)[0] / df['PM2.5'])
#    print(mass_closure_fraction)
    #for i in range(0, len(mass_closure_fraction)):
    #    print(mass_closure_fraction[0])
         #print('Filtro ' + str(i + 1) + ' = ' + str(mass_closure_fraction[i]))


        
for i in methods:
    print('\n\n'+ str(i))
    mass = mass_closure(data_df=df, equation=i)
    print('Total explained mass = ', ((mass[0] / df['PM2.5']).mean()* 100).round(1), '%')
    for key in mass[1].keys():
        print(key, '=', ((mass[1][key] / df['PM2.5']).mean() * 100).round(1), '%')



mass = mass_closure(data_df=df, equation='Hand_2011')[1]
labels = []
sizes = []

for x, y in mass.items():
    labels.append(x)
    sizes.append(y.mean())


print(labels)
labels_sty = ['II', 'OM', 'EC', 'GM', 'SSA', 'S/E']
# Plot
plt.figure(figsize=(4.5,5))
plt.pie(sizes, labels=labels_sty)

plt.axis('equal')
plt.savefig('imagenes/torta_media_hand.pdf')
plt.show()
# closure, inorganic_ions, organic_mass, elemental_C,
#geological_minerals, salt, trace_elements, others, unexplained =




# +


methods = ['Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011', 'Hand_2011_mod']

labels = ['Solomon 1989', 'Chow 1994', 'Malm 1994', 'Chow 1996', 'Andrews 2000',
         'Malm 2000', 'Maenhaut 2002', 'DeBell 2006', 'Hand 2011', 'Hand 2011 mod']

for m in methods:
    escapan = 0
    escapan_defecto = 0
    escapan_exceso = 0
    cumplen = 0
    son_nan = 0
    criterio = 0.2
    for i in range(0,120-38):
        value = ((mass_closure(data_df=df, equation=m)[0]) / df['PM2.5']).iloc[i]
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


style='.-'
figure, ax = plt.subplots()
hand = ((mass_closure(data_df=df, equation='Hand_2011')[0] / df['PM2.5'])).plot(style=style, label='Hand 2011')
hand_mod = (mass_closure(data_df=df, equation='Hand_2011_mod')[0] / df['PM2.5']).plot(style=style, color='r', label='Hand 2011 mod')
            
#((mass_closure(data_df=df, equation='DeBell_2006')[0] / df['PM2.5']).plot(style='o-', color='g', label='DeBell'))
#((mass_closure(data_df=df, equation='Malm_2000')[0] / df['PM2.5']).plot(style='o-', color='k', label='Malm'))
plt.legend()

plt.hlines(.8, xmin=df.index[0], xmax=df.index[-1], color='k', linestyles=':')
plt.hlines(1.2, xmin=df.index[0], xmax=df.index[-1], color='k', linestyles=':')
#plt.xticks(np.arange(30,120, 10))
#plt.grid()
plt.savefig('imagenes/reconstruccion_masica_hand2011.pdf')
plt.show()




# +
methods = ['Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011', 'Hand_2011_mod']

labels = ['Solomon 1989', 'Chow 1994', 'Malm 1994', 'Chow 1996', 'Andrews 2000',
         'Malm 2000', 'Maenhaut 2002', 'DeBell 2006', 'Hand 2011', 'Hand 2011 mod']

figs, axs = plt.subplots(2, 5, figsize=(16,5))
for i in range(len(methods) - 1):
    style='.-'
    xax = i%2
    yax = int(i/2)
    axs[xax, yax].plot((mass_closure(data_df=df, equation=methods[i])[0] / df['PM2.5']), label=labels[i])
    axs[xax, yax].set_title(labels[i])

    axs[xax, yax].hlines(.8, xmin=df.index[0], xmax=df.index[-1], color='k', linestyles=':')
    axs[xax, yax].hlines(1.2, xmin=df.index[0], xmax=df.index[-1], color='k', linestyles=':')

for ax in axs.flat:
    ax.set(xlabel='Date', ylabel='Fraction explained')
    ax.tick_params(labelrotation=45)

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
#plt.setp(axs.xaxis.get_majorticklabels(), rotation=45)


# -

value = ((mass_closure(data_df=df, equation='Hand_2011')[0]) / df['PM2.5'])[41]
print(np.isnan(value))

# + tags=[]
# Tortas
methods = ['Solomon_1989', 'Chow_1994', 'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006', 'Hand_2011']

#sizes = (mass_closure(data_df=df.iloc[1], equation=k)[1:])
for k in methods:
    print(k)
    for j in range(len(df)):
        print(df.index[j])
        sizes = (mass_closure(data_df=df.iloc[j], equation=k)[1])
        labels = list(sizes.keys())
        values = list(sizes.values())
        try:
            fig1, ax1 = plt.subplots()
            ax1.pie(values, labels=labels)
            ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
            plt.savefig('imagenes/tortas/torta_filtro_' + 
                        str(df.index[j]) + '_metodo_' + k + '.eps')
            plt.show()
        except:
            print('no se pudo hacer', df.index[j])

# + tags=[]
# Correlaciones de todo con todo
graph_all_corr(df, 'correlaciones/')
# -
(df['C Elemental'] / df['C Orgánico']).plot(style='.-', label="EC/OC")
plt.legend()
plt.show()

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

# + tags=[]
# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(df['NO3'], df['Cl'], df['Na'], 'o')
# ax.set_xlabel(r'NO$_3^-$')
# ax.set_ylabel(r'Cl$^-$')
# ax.set_zlabel(r'Na$^+$')
# plt.show()
# 
# plt.figure(figsize=(10,8))
# plt.scatter(df['NO3'], df['Cl'], c=df['Na'], cmap='viridis')
# plt.xlabel('NO3')
# plt.xlim(0, 1.5)
# plt.ylim(0, 0.4)
# plt.ylabel('Cl')
# plt.colorbar()
# plt.show()
# 
# 
# plt.figure(figsize=(10,8))
# plt.scatter(df['Na'], df['NO3'], c=np.log10(df['Cl']), cmap='viridis')
# plt.xlabel('Na')
# #plt.xlim(0, 1.5)
# plt.ylim(0, 1.5)
# plt.ylabel('NO3')
# plt.colorbar()
# plt.show()
# 
valores_cloro = np.arange(0.4, 0.0, -0.05)

fig, axs = plt.subplots(4, 2, figsize=(10, 20))
iterador = 0
for i in valores_cloro:
    limite = i.round(3)
    cond1 = df['Cl'] >= limite
    cond2 = df['NO3'] <= 1.5
    ax = axs[iterador // 2, iterador % 2]
    ax.plot(df['Na'].where(cond1).where(cond2), df['NO3'].where(cond1).where(cond2), 'o', label=r'Cl $\geq$ ' + str(limite))
    ax.set_ylim(0, 1.5)
    ax.set_xlim(0, 0.7)
    ax.set_xlabel('Na$^+$')
    ax.set_ylabel('NO$_3^-$')
    ax.legend()
    iterador += 1
plt.show()
