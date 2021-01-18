# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.7.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
from pandas.api.types import is_numeric_dtype
import xarray as xr

# %run ./funciones.ipynb

# +
df_raw = pd.read_csv('20210105_matriz_medidas_valores_apreciables.csv', index_col=0, skiprows=[1,2])[0:94]
#df['PM2.5 (ug/m3)'] = df['PM2.5 (mg/m3)']*1000 * (df['PM2.5 (mg/m3)'] < 1) * (df['PM2.5 (mg/m3)'] > 0)

############# Remove outliers ############
df = df_raw.mask(df_raw.sub(df_raw.mean()).div(df_raw.std()).abs().gt(3))

df_tecnicas = pd.read_csv('ArcalMetalesAnalisis.csv', index_col=0, skiprows=[1,2])
# -


#print(df_tecnicas['Ti'])
# %run ./funciones.ipynb
lista_elementos = ['Ti', 'Mn', 'Ni', 'Zn', 'Pb']
for i in lista_elementos:
    ax = comparacion_tecnicas(df_tecnicas, i)
    plt.show()

# +
#display(df.corr())
df.corr().to_csv('correlacion_todo.csv')

df.plot(y=['C Orgánico', 'C Elemental', 'C Total', 'PM2.5'], ylim=(0,20))
df.hist('C Orgánico', bins=20)
df.hist('C Elemental', bins=20)
df.hist('C Total', bins=20)
plt.show()
log_PM25 = np.log(df['PM2.5'].where(df['PM2.5']>0.1))
log_PM25.hist(bins=20)
plt.title('ln(PM2.5)')
plt.show()

### This function tests the null hypothesis that a sample comes from a normal distribution.
### It is based onD’Agostino and Pearson’s [1], [2] test that combines skew and kurtosis 
### to produce an omnibus test of normality.
print('Test de normalidad ln(PM2.5), p-valor:')
print(stats.normaltest(log_PM25, axis=0, nan_policy='omit').pvalue)
#df.plot(y=['C Elemental (ug/cm2)'], xlim=(0,36), ylim=(0,0.75))

df.hist('SO4', bins=20)
plt.show()
log_SO4 = np.log(df['SO4'])
log_SO4.hist(bins=20)
plt.title('ln(SO4)')
plt.show()

df.hist('Pb', bins=20)
plt.show()
log_SO4 = np.log(df['Pb'])
log_SO4.hist(bins=20)
plt.title('ln(Pb)')
plt.show()

df.hist('Ca', bins=20)
plt.show()
log_SO4 = np.log(df['Ca'])
log_SO4.hist(bins=20)
plt.title('ln(Ca)')
plt.show()
# -

not_lognormal = []
for i in df:
    if not df[i].isna().all() and i != 'Se.1': 
        p_is_lognormal = test_lognormality(df[i])
        print(i, str(p_is_lognormal))
        if p_is_lognormal < 0.01:
            not_lognormal.append(i)
print('\n\n')
print(not_lognormal)

# +
#for i in df:
#    if is_numeric_dtype(df[i]):
#        #plt.figure(figsize=(20,10))
#        df.plot(y=i, label=i, xlim=(0, 96))
#        plt.show()
# -

#df['C Orgánico'].plot(label='Orgánico')
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
#df['C Orgánico'].plot(label='Orgánico')
fraccion_org = (df['C Orgánico']*1.4) / df['PM2.5']
fraccion_elemental = df['C Elemental'] / df['PM2.5']
(fraccion_org).plot(label='C Orgánico/PM2.5')
plt.legend()
plt.show()
print(fraccion_org.mean())

#print(fraccion_elemental.mean())
# -

fraccion_NO3 = df['NO3'] / df['PM2.5']
(fraccion_NO3).plot(label='NO3/PM2.5')
plt.legend()
plt.show()
print(fraccion_NO3.mean())

# +
for i in df:
    if 'KED' not in i:
        print(i)
        print(df[i].mean()/df['PM2.5'].mean())

masa_explicada = 0
for i in df:
    if '.1' not in i:
        masa_explicada += df[i].mean()
masa_explicada = (masa_explicada - df['PM2.5'].mean() - df['C Total'].mean())/(df['PM2.5'].mean())
print('masa explicada = ', str(masa_explicada))
# -

# (df['C Total'] *1.4 * 2.5).plot(label='Orgánico')
# (df['Na'] * 10).plot(label='Na')
# (df['S']).plot(label='S', xlim=(0,40))
# (df['NO3']).plot(label='NO$_3$')
# (df['NH4']).plot(label='NH$_4$')
# plt.legend()
# plt.show()
#
# df['S'].plot(label='S', xlim=(0,40))
# df['SO4'].plot(label='SO$_4$')
# plt.legend()
# plt.show()
#
# df['Ni.1'].plot(label='Ni.1')
# df['Ni'].plot(label='Ni')
# plt.legend()
# plt.show()
#
# df['Pb'].plot(label='Pb')
# df['Sb'].plot(label='Sb')
# plt.legend()
# plt.show()



# +
#especies_corr_int = [['NO3', 'S', 'PM2.5'],
#                          ['S', 'Ca', 'Zn', 'PM2.5', 'NH4'], ['Na', 'Fe'],
#                          ['C Orgánico', 'Sb'], ['Ca', 'Sr'], ['Cr', 'Ni', 'Ni.1'],
#                          ['Mn', 'Mn.1'], ['Pb', 'Zn'], ['Pb.1', 'Sb'],
#                          ['Mn', 'Ti']]
#for i in especies_corr_int:
#    for j in i:
#        if j == 'PM2.5':
#            (df[j]*0.1).plot(label=j + ' * 0.1', xlim=(0,40))
#        else:
#            df[j].plot(label=j, xlim=(0,40))
#    plt.legend()
#    plt.show()
# -

#df.corr().where(df.corr() > 0.8).to_csv('correlacion_grande.csv')
corr_matrix(df, minimum=0.8).to_csv('correlacion_grande.csv')

# +
# Correlaciones de todo con todo
# graph_all_corr(df, correlaciones)
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
    mass_closure_fraction = (mass_closure(data_df=df[0:40], equation=k) / df['PM2.5'][0:40])
    print(mass_closure_fraction)
    #for i in range(0, len(mass_closure_fraction)):
    #    print(mass_closure_fraction[0])
         #print('Filtro ' + str(i + 1) + ' = ' + str(mass_closure_fraction[i]))


for i in methods:
    print('\n\n',i)
    print(str((mass_closure(data_df=df[0:40], equation=i) / df['PM2.5'][0:40]).mean()))
        
        
