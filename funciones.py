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

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
from pandas.api.types import is_numeric_dtype
import xarray as xr


def corr_matrix(df, minimum=-1, maximum=1):
    '''Calcula matriz de correlación y devuelve los valores entre las corr
    máximas y mínimas'''
    return df.corr().where(df.corr() >= minimum).where(df.corr() <= maximum)


def test_lognormality(data_array, treshold=0):
    '''testea si es lognormal usando un método basado en D'Agostino y Pearson,
    que tiene en cuenta la distosión y kurtosis para producir el test'''
    try:
        pvalor = stats.normaltest(np.log(data_array.where(
                                  data_array > treshold)), axis=0, nan_policy='omit').pvalue
    except:
        pvalor = np.nan
    return pvalor


def graph_all_corr(df, output_dir='./'):
    '''Grafica correlación de todo con todo y la imprime'''
    keys = df.keys()
    for i in range(0, len(keys)):
        if (not df[keys[i]].isna().all()) and is_numeric_dtype(df[keys[i]]):
            for j in range(i + 1, len(keys)):
                print('i == ' + keys[i] + ', j == ' + keys[j])
                if (not df[keys[j]].isna().all()) and is_numeric_dtype(df[keys[j]]):
                    try:
                        idx = np.isfinite(df[keys[i]]) & np.isfinite(df[keys[j]])
                        coef = np.polyfit(df[keys[i]][idx],df[keys[j]][idx],1)
                        print(str(coef[0]), 'x +', str(coef[1]))
                        poly1d_fn = np.poly1d(coef)
                    except:
                        print('No se pudo hacer i ==' + keys[i] + ', j == ' + keys[j])
                    else:
                        correlation_fig = plt.plot(df[keys[i]],df[keys[j]], 'yo',
                                                   df[keys[i]], poly1d_fn(df[keys[i]]), 'k')
                        plt.legend(iter(correlation_fig), ('Values', 'Linear fit'))
                        #plt.plot(df[i], df[i], 'r--')
                        plt.xlabel(keys[i])
                        plt.ylabel(keys[j])
                        plt.savefig(output_dir + 'corr_' + keys[i] + '_' + keys[j] + '.png')
                        plt.show()


def comparacion_tecnicas(df, elemento):
    try:
        idx = np.isfinite(df[elemento]) & np.isfinite(df[elemento + '.1'])
        coef = np.polyfit(df[elemento][idx], df[elemento + '.1'][idx],1)
        coeficients = str(coef[0]), 'x + ', str(coef[1])
        poly1d_fn = np.poly1d(coef)
    except:
        print('No se pudo hacer i == ' + elemento)
    else:
        fig, ax = plt.subplots()
        correlation_fig = plt.plot(df[elemento], df[elemento + '.1'],
                                   'yo', df[elemento], poly1d_fn(df[elemento]), 'k', axes=ax)
        ax.text(0.1, 0.9,
                str(coef[0].round(2)) +  'x + ' + str(coef[1].round(4)), transform=ax.transAxes)
        plt.xlabel(elemento + ' (TXRF)')
        plt.ylabel(elemento + ' (ICP-MS)')
        return(ax)


# +
# df.keys() = ['Cl', 'NO3', 'SO4', 'Na', 'NH4', 'C Orgánico', 'C Elemental', 'C Total',
       # 'S', 'Cl.1', 'K', 'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Ni', 'Cu', 'Zn',
       # 'As', 'Se', 'Sr', 'Pb', 'Hg', 'Pb.1', 'Zn.1', 'Mn.1', 'Mo', 'Ni.1',
       # 'Ti.1', 'Sb', 'PM2.5']
    


def mass_closure(data_df, equation='Chow_1996'):
    closure = 0
    data_df2 = data_df.fillna(0)
    data_df2['Si'] = 2.4729 * data_df2['Al']
    
    if equation == 'Macias_1981':
        inorganic_ions = data_df2['(NH4)2SO4'] + data_df2['NH4SO3']
        organic_mass = 1.5 * data_df2['C Orgánico']
        elemental_C = data_df2['C Elemental']
        geological_minerals = (1.89 * data_df2['Al'] + 2.14 * data_df2['Si'] +
                               1.4 * data_df2['Ca'] + 1.2 * data_df2['K'] +
                               1.43 * data_df2['Fe'])
        trace_elements = (1.25 * data_df2['Cu'] + 1.24 * data_df2['Zn'] +
                          1.08 * data_df2['Pb'])
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                       'trace_elements': trace_elements}
      
    if equation == 'Solomon_1989':
        inorganic_ions = data_df2['SO4'] + data_df2['NO3'] + data_df2['NH4']
        organic_mass = 1.4 * data_df2['C Orgánico']
        elemental_C = data_df2['C Elemental']
        geological_minerals = (1.89 * data_df2['Al'] + 2.14 * data_df2['Si'] +
                               1.4 * data_df2['Ca'] + 1.43 * data_df2['Fe'])
        # Pide que sean por XRF salvo Na+ y Mg++, que deberían ser por AAS. No es el caso en ninguno
        trace_elements = (data_df2['Cl'] + data_df2['Na'] + data_df2['K'] +
                          data_df2['Ti'] + data_df2['V'] + data_df2['Cr'] +
                          data_df2['Mn'] + data_df2['Ni'] + data_df2['Cu'] +
                          data_df2['Zn'] + data_df2['As'] + data_df2['Se'] +
                          data_df2['Sr'] + data_df2['Pb'] + data_df2['Hg'] +
                          data_df2['Sb'])
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'trace_elements': trace_elements}
        
    if equation == 'Chow_1994':
        inorganic_ions = data_df2['SO4'] + data_df2['NO3'] + data_df2['NH4']
        organic_mass = 1.4 * data_df2['C Orgánico']
        elemental_C = data_df2['C Elemental']
        geological_minerals = (1.89 * data_df2['Al'] + 2.14 * data_df2['Si'] +
                               1.4 * data_df2['Ca'] + 1.43 * data_df2['Fe'])
        trace_elements = (data_df2['Cl'] + data_df2['Na'] + data_df2['K'] +
                          data_df2['Ti'] + data_df2['V'] + data_df2['Cr'] +
                          data_df2['Mn'] + data_df2['Ni'] + data_df2['Cu'] +
                          data_df2['Zn'] + data_df2['As'] + data_df2['Se'] +
                          data_df2['Sr'] + data_df2['Pb'] + data_df2['Hg'] +
                          data_df2['Sb'])
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'trace_elements': trace_elements}
        
    if equation == 'Malm_1994':
        inorganic_ions = 4.125 * data_df2['S']
        organic_mass = 1.4 * data_df2['C Orgánico']
        elemental_C = data_df2['C Elemental']
        geological_minerals = (2.2 * data_df2['Al'] + 2.49 * data_df2['Si'] +
                               1.63 * data_df2['Ca'] + 1.94 * data_df2['Ti'] +
                               2.42 * data_df2['Fe'])
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals}
        
    if equation == 'Chow_1996':
        inorganic_ions = data_df2['SO4'] + data_df2['NO3'] + data_df2['NH4']
        organic_mass = 1.4 * data_df2['C Orgánico']
        elemental_C = data_df2['C Elemental']
        geological_minerals = (1.89 * data_df2['Al'] + 2.14 * data_df2['Si'] +
                               1.4 * data_df2['Ca'] + 1.43 * data_df2['Fe'])
        salt = data_df2['Na'] + data_df2['Cl']
        trace_elements = (data_df2['Ti'] + data_df2['V'] + data_df2['Cr'] +
                          data_df2['Mn'] + data_df2['Ni'] + data_df2['Cu'] +
                          data_df2['Zn'] + data_df2['As'] + data_df2['Se'] +
                          data_df2['Sr'] + data_df2['Pb'] + data_df2['Hg'] +
                          data_df2['Sb'])
       
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'salt': salt, 'trace_elements': trace_elements}
        
    if equation == 'Andrews_2000':
        inorganic_ions = data_df2['SO4'] + data_df2['NO3'] + data_df2['NH4']
        organic_mass = 1.4 * data_df2['C Orgánico']
        elemental_C = data_df2['C Elemental']
        geological_minerals = (1.89 * data_df2['Al'] + 2.14 * data_df2['Si'] +
                               1.4 * data_df2['Ca'] + 1.43 * data_df2['Fe'] +
                               1.67 * data_df2['Ti'])
        trace_elements = (data_df2['Cl'] + data_df2['Na'] + data_df2['K'] +
                          data_df2['V'] + data_df2['Cr'] +
                          data_df2['Mn'] + data_df2['Ni'] + data_df2['Cu'] +
                          data_df2['Zn'] + data_df2['As'] + data_df2['Se'] +
                          data_df2['Sr'] + data_df2['Pb'] + data_df2['Hg'] +
                          data_df2['Sb'])
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'trace_elements': trace_elements}
        
    if equation == 'Malm_2000':
        inorganic_ions = 1.125 * data_df2['S'] + 1.29 * data_df2['NO3']
        organic_mass = 1.4 * data_df2['C Orgánico']
        elemental_C = data_df2['C Elemental']
        geological_minerals = (2.2 * data_df2['Al'] + 2.49 * data_df2['Si'] +
                               1.63 * data_df2['Ca'] + 1.94 * data_df2['Ti'] +
                               2.42 * data_df2['Fe'])
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals}
    
    if equation == 'Maenhaut_2002':
        inorganic_ions = data_df2['SO4'] + data_df2['NO3'] + data_df2['NH4']
        organic_mass = 1.4 * data_df2['C Orgánico']
        elemental_C = data_df2['C Elemental']
        geological_minerals = (2.2 * data_df2['Al'] + 2.49 * data_df2['Si'] +
                               1.63 * data_df2['Ca'] + 1.94 * data_df2['Ti'] +
                               2.42 * data_df2['Fe'])
        salt = data_df2['Cl'] + 1.4486 * data_df2['Na']
        trace_elements =(data_df2['K'] +  data_df2['V'] + data_df2['Cr'] +
                          data_df2['Mn'] + data_df2['Ni'] + data_df2['Cu'] +
                          data_df2['Zn'] + data_df2['As'] + data_df2['Se'] +
                          data_df2['Sr'] + data_df2['Pb'] + data_df2['Hg'] +
                          data_df2['Sb'])  # CHEQUEAR
        others = data_df['K'] - 0.6 * data_df['Fe']
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
              'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
              'salt': salt, 'trace_elements': trace_elements, 'others': others}
        
    if equation == 'DeBell_2006':
        inorganic_ions = 4.125 * data_df2['S'] + 1.29 * data_df2['NO3']
        organic_mass = 1.8 * data_df2['C Orgánico']
        elemental_C = data_df2['C Elemental']
        geological_minerals = (2.2 * data_df2['Al'] + 2.49 * data_df2['Si'] +
                               1.63 * data_df2['Ca'] + 1.94 * data_df2['Ti'] +
                               2.42 * data_df2['Fe'])
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals}
        
    if equation == 'Hand_2011':
        inorganic_ions = 1.375 * data_df2['SO4'] + 1.29 * data_df2['NO3']
        organic_mass = 1.8 * data_df2['C Orgánico']
        elemental_C = data_df2['C Elemental']
        geological_minerals = (2.2 * data_df2['Al'] + 2.49 * data_df2['Si'] +
                               1.63 * data_df2['Ca'] + 1.94 * data_df2['Ti'] +
                               2.42 * data_df2['Fe'])
        salt = 1.8 * data_df2['Cl']
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'salt': salt}
        
    if equation == 'Hand_2011_mod':
        inorganic_ions = 1.375 * data_df2['SO4'] + 1.29 * data_df2['NO3']
        # organic_mass = 1.8 * data_df2['C Orgánico']
        organic_mass = 2.08 * data_df2['C Orgánico']
        elemental_C = data_df2['C Elemental']
        geological_minerals = (2.2 * data_df2['Al'] + 2.49 * data_df2['Si'] +
                               1.63 * data_df2['Ca'] + 1.94 * data_df2['Ti'] +
                               2.42 * data_df2['Fe'])
        salt = 1.8 * data_df2['Cl']

        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
              'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
              'salt': salt}
        
    if equation == 'Simon_2011':
        inorganic_ions = data_df2['(NH4)2SO4'] + data_df2['NH4NO3']
        organic_mass = 1.8 * data_df2['C Orgánico']
        elemental_C = data_df2['C Elemental']
        geological_minerals = (3.48 * data_df2['Si'] + 1.63 * data_df2['Ca'] +
                               2.42 * data_df2['Fe'] + 1.94 * data_df2['Ti'])
        salt = 1.8 * data_df2['Cl']
        others = 1.2 * (data_df['K'] - 0.6 * data_df['Fe'])
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
              'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
              'salt': salt, 'others': others}
        
    closure = sum(categories.values())
    categories['unexplained'] = data_df2['PM2.5'] - closure
    return (closure, categories)



# -

class Filtro:
    '''Clase que contiene los 
    datos relevantes de los filtros'''
    def __init__(self, da):
        self.mass = da['PM2.5']
        self.Cl = da['Cl']
        self.NO3 = da['NO3']
        self.Na = da['Na']
        self.NH4 = da['NH4']
        self.SO4 = da['SO4']
    
    def Solomon_1989(self, value='inorganic_ions'):
        if value == 'inorganic_ions':
            return self.SO4 + self.NO3 + self.NH4
        
    def Chow_1994(self, value='inorganic_ions'):
        if value == 'inorganic_ions':
            return self.SO4 + self.NO3 + self.NH4
        
    def Chow_1996(self, value='inorganic_ions'):
        if value == 'inorganic_ions':
            return self.SO4 + self.NO3 + self.NH4
        
    def Andews_2000(self, value='inorganic_ions'):
        if value == 'inorganic_ions':
            return self.SO4 + self.NO3 + self.NH4
        
    def Maenhaut_2002(self, value='inorganic_ions'):
        if value == 'inorganic_ions':
            return self.SO4 + self.NO3 + self.NH4
        
    def Hand_2011(self, value='inorganic_ions'):
        if value == 'inorganic_ions':
            return 1.375 * self.SO4 + 1.29 * self.NO3
