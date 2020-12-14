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
df_raw = pd.read_csv('04122020_matriz_medidas_valores_apreciables.csv', index_col=0, skiprows=[1,2])
#df['PM2.5 (ug/m3)'] = df['PM2.5 (mg/m3)']*1000 * (df['PM2.5 (mg/m3)'] < 1) * (df['PM2.5 (mg/m3)'] > 0)

############# Remove outliers ############
df = df_raw.mask(df_raw.sub(df_raw.mean()).div(df_raw.std()).abs().gt(3))
display(df)


# +
display(df.corr())
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
# -

not_lognormal = []
for i in df:
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
    print(i)
    print(df[i].mean()/df['PM2.5'].mean())

masa_explicada = 0
for i in df:
    if '.1' not in i:
        masa_explicada += df[i].mean()
masa_explicada = (masa_explicada - df['PM2.5'].mean() - df['C Total'].mean())/(df['PM2.5'].mean())
print('masa explicada = ', str(masa_explicada))

# +
#(df['C Total'] *1.4 * 2.5).plot(label='Orgánico')
#(df['Na'] * 10).plot(label='Na')
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

especies_corr_int = [['NO3', 'S', 'PM2.5'],
                          ['S', 'Ca', 'Zn', 'PM2.5', 'NH4'], ['Na', 'Fe'],
                          ['C Orgánico', 'Sb'], ['Ca', 'Sr'], ['Cr', 'Ni', 'Ni.1'],
                          ['Mn', 'Mn.1'], ['Pb', 'Zn'], ['Pb.1', 'Sb'],
                          ['Mn', 'Ti']]
for i in especies_corr_int:
    for j in i:
        if j == 'PM2.5':
            (df[j]*0.1).plot(label=j + ' * 0.1', xlim=(0,40))
        else:
            df[j].plot(label=j, xlim=(0,40))
    plt.legend()
    plt.show()

#df.corr().where(df.corr() > 0.8).to_csv('correlacion_grande.csv')
corr_matrix(df, minimum=0.8).to_csv('correlacion_grande.csv')

graph_all_corr(df)

plt.plot((df['S'] - df['SO4'] * 0.33), df['V'], 'yo')
plt.show()
(df['V'] * 100).plot()
(df['S'] - df['SO4'] * 0.33).plot()


# +
def mass_closure(data_df, method='Chow_1996'):
    mass_closure = 0
    data_df2 = data_df.fillna(0)
    if method == 'Chow_1996':
        ##### FALTA EXCEPTUAR ALUMINIO, Silicio
        data_df2['Al'] = 0 * data_df2['Ca']
        data_df2['Si'] = 0 * data_df2['Ca']
        ########
        inorganic_ions = data_df2['SO4'] + data_df2['NO3'] + data_df2['NH4']
        organic_matter = data_df2['C Orgánico'] * 1.4
        elemental_C = data_df2['C Elemental']
        geological_minerals = (1.89 * data_df2['Al'] + 2.14 * data_df2['Si'] +
                               1.4 * data_df2['Ca'] + 1.43 * data_df2['Fe'])
        salt = data_df2['Na'] + data_df2['Cl']
        trace_elements = data_df2['V'] #Terminar
        others = data_df2['Na'] * 0
    
    
    mass_closure = (inorganic_ions + organic_matter + elemental_C + geological_minerals +
                    salt + trace_elements + others)
    return mass_closure

for i in (mass_closure(data_df=df) / df['PM2.5']):
    print(i)
    
print('\n\n\n' + str((mass_closure(data_df=df) / df['PM2.5']).mean()))
        
        ######### Metodo de Malm 1994
        #geological_minerals = (2.20 * data_df['Al'] + 2.49 * data_df['Si'] +
        #                       1.63 * data_df['Ca'] + 1.94 * data_df['Ti'] +
        #                       2.42 * data_df['Fe'])
        
