#!/usr/bin/env python
# coding: utf-8

# In[ ]:
 

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
from pandas.api.types import is_numeric_dtype
import xarray as xr
import statsmodels.api as sm
from scipy.stats import linregress



def corr_matrix(df, minimum=-1, maximum=1, method='pearson'):
    '''Calcula matriz de correlación y devuelve los valores entre las corr
    máximas y mínimas'''
    correlation = df.corr(method)
    return correlation.where(correlation >= minimum).where(correlation <= maximum)


def mass_reconstruction(conc_matrix, unc_matrix, equation='Hand_2011'):
    
    """
    Reconstructs the mass using methodologies as described in Chow 2015. It requires a concentration
    and an uncertanty matrices.
    """

    closure = 0
    concentration_matrix = conc_matrix.copy()
    uncertainty_matrix = unc_matrix.copy()
#    data_df2 = data_df.fillna(0)
    if 'Si' not in concentration_matrix:
        concentration_matrix['Si'] = 2.4729 * concentration_matrix['Al']
        uncertainty_matrix['Si'] = 2.4729 * uncertainty_matrix['Al'] * 2 # Ese * 2 es solo para agrandar la incert Si
    concentration_matrix['S'] = concentration_matrix['SO4']/3 
    uncertainty_matrix['S'] = uncertainty_matrix['SO4']/3 
    concentration_matrix['Sr'] = 0
    uncertainty_matrix['Sr'] = 0
    concentration_matrix['Hg'] = 0
    uncertainty_matrix['Hg'] = 0
    
    if equation == 'Macias_1981':

                
        if "(NH4)2SO4" not in concentration_matrix:
            # Assume all SO4 is (NH4)2SO4
            concentration_matrix["(NH4)2SO4"] = (132.14 / 96.06 ) * concentration_matrix["SO4"]
            uncertainty_matrix["(NH4)2SO4"] = (132.14 / 96.06 ) * uncertainty_matrix["SO4"]
        if "NH4NO3" not in concentration_matrix:
            # Assume all NO3 is NH4NO3
            concentration_matrix["NH4NO3"] =  (80.043 / 62.004 ) * concentration_matrix["NO3"]
            uncertainty_matrix["NH4NO3"] = (80.043 / 62.004 ) * uncertainty_matrix["NO3"]
        
        
        inorganic_ions = concentration_matrix['(NH4)2SO4'] + concentration_matrix['NH4NO3']
        uinorganic_ions = np.linalg.norm( [uncertainty_matrix['(NH4)2SO4'], uncertainty_matrix['NH4NO3'] ], axis=0)
        
        organic_mass = 1.5 * concentration_matrix['OC']
        uorganic_mass = 1.5 * uncertainty_matrix['OC']
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (1.89 * concentration_matrix['Al'] + 2.14 * concentration_matrix['Si'] +
                               1.4 * concentration_matrix['Ca'] + 1.2 * concentration_matrix['K'] +
                               1.43 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                                1.4 * uncertainty_matrix['Ca'], 1.2 * uncertainty_matrix['K'],
                                                1.43 * uncertainty_matrix['Fe'] ], axis=0)
        
        trace_elements = (1.25 * concentration_matrix['Cu'] + 1.24 * concentration_matrix['Zn'] +
                          1.08 * concentration_matrix['Pb'])
        utrace_elements = np.linalg.norm( [1.25 * uncertainty_matrix['Cu'], 1.24 * uncertainty_matrix['Zn'],
                                              1.08 * uncertainty_matrix['Pb'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                       'trace_elements': trace_elements}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                       'utrace_elements': utrace_elements}
        
        uclosure = np.linalg.norm( [ uncertainty_matrix['(NH4)2SO4'], uncertainty_matrix['NH4NO3'],
                                        1.5 * uncertainty_matrix['OC'],
                                        uncertainty_matrix['EC'],
                                        1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                        1.4 * uncertainty_matrix['Ca'], 1.2 * uncertainty_matrix['K'],
                                        1.43 * uncertainty_matrix['Fe'],
                                        1.25 * uncertainty_matrix['Cu'], 1.24 * uncertainty_matrix['Zn'],
                                        1.08 * uncertainty_matrix['Pb'] ], axis=0)
                               
      
    if equation == 'Solomon_1989':
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm([ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'] ], axis=0)
        
        organic_mass = 1.4 * concentration_matrix['OC']
        uorganic_mass = 1.4 * uncertainty_matrix['OC']
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (1.89 * concentration_matrix['Al'] + 2.14 * concentration_matrix['Si'] +
                               1.4 * concentration_matrix['Ca'] + 1.43 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [1.89 * uncertainty_matrix['Al'] + 2.14 * uncertainty_matrix['Si'] +
                               1.4 * uncertainty_matrix['Ca'] + 1.43 * uncertainty_matrix['Fe'] ], axis=0)
        
        # Pide que sean por XRF salvo Na+ y Mg++, que deberían ser por AAS. No es el caso en ninguno
        trace_elements = (concentration_matrix['Cl'] + concentration_matrix['Na sol'] + concentration_matrix['Na no sol'] +
                          concentration_matrix['K'] +
                          concentration_matrix['Ti'] + concentration_matrix['V'] + concentration_matrix['Cr'] +
                          concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                          concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                          concentration_matrix['Sr'] + concentration_matrix['Pb'] + concentration_matrix['Hg'] +
                          concentration_matrix['Sb'])
        utrace_elements = np.linalg.norm( [uncertainty_matrix['Cl'] + uncertainty_matrix['Na sol'] + uncertainty_matrix['Na no sol'] +
                                           uncertainty_matrix['K'] +
                          uncertainty_matrix['Ti'] + uncertainty_matrix['V'] + uncertainty_matrix['Cr'] +
                          uncertainty_matrix['Mn'] + uncertainty_matrix['Ni'] + uncertainty_matrix['Cu'] +
                          uncertainty_matrix['Zn'] + uncertainty_matrix['As'] + uncertainty_matrix['Se'] +
                          uncertainty_matrix['Sr'] + uncertainty_matrix['Pb'] + uncertainty_matrix['Hg'] +
                          uncertainty_matrix['Sb'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'trace_elements': trace_elements}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                      'utrace_elements': utrace_elements}
        
        uclosure = np.linalg.norm( [uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'],
                                    1.4 * uncertainty_matrix['OC'],
                                    uncertainty_matrix['EC'],
                                    1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                    1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'],
                                    uncertainty_matrix['Cl'] + uncertainty_matrix['Na sol'] + uncertainty_matrix['Na no sol'] +
                                    uncertainty_matrix['K'] +
                                    uncertainty_matrix['Ti'] + uncertainty_matrix['V'] + uncertainty_matrix['Cr'] +
                                    uncertainty_matrix['Mn'] + uncertainty_matrix['Ni'] + uncertainty_matrix['Cu'] +
                                    uncertainty_matrix['Zn'] + uncertainty_matrix['As'] + uncertainty_matrix['Se'] +
                                    uncertainty_matrix['Sr'] + uncertainty_matrix['Pb'] + uncertainty_matrix['Hg'] +
                                    uncertainty_matrix['Sb'] ], axis=0)
                                    
        
    if equation == 'Chow_1994':
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'] ], axis=0)
        
        organic_mass = 1.4 * concentration_matrix['OC']
        uorganic_mass = 1.4 * uncertainty_matrix['OC']
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (1.89 * concentration_matrix['Al'] + 2.14 * concentration_matrix['Si'] +
                               1.4 * concentration_matrix['Ca'] + 1.43 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'] ], axis=0)
        
        trace_elements = (concentration_matrix['Cl'] + concentration_matrix['Na sol'] + concentration_matrix['Na no sol'] +
                          concentration_matrix['K'] +
                          concentration_matrix['Ti'] + concentration_matrix['V'] + concentration_matrix['Cr'] +
                          concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                          concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                          concentration_matrix['Sr'] + concentration_matrix['Pb'] + concentration_matrix['Hg'] +
                          concentration_matrix['Sb'])
        utrace_elements = np.linalg.norm( [ uncertainty_matrix['Cl'], uncertainty_matrix['Na sol'],  uncertainty_matrix['Na no sol'],
                                            uncertainty_matrix['K'],
                                            uncertainty_matrix['Ti'], uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                                            uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                            uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                            uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                                            uncertainty_matrix['Sb'] ], axis=0)     
        
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'trace_elements': trace_elements}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                       'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                       'utrace_elements': utrace_elements}
        
        uclosure = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'],
                                    1.4 * uncertainty_matrix['OC'],
                                    uncertainty_matrix['EC'],
                                    1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                    1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'],
                                    uncertainty_matrix['Cl'], uncertainty_matrix['Na sol'], uncertainty_matrix['Na no sol'],
                                    uncertainty_matrix['K'],
                                    uncertainty_matrix['Ti'], uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                                    uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                    uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                    uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                                    uncertainty_matrix['Sb'] ], axis=0)
                                    
        
    if equation == 'Malm_1994':
        inorganic_ions = 4.125 * concentration_matrix['S']
        uinorganic_ions = 4.125 * uncertainty_matrix['S']
        
        organic_mass = 1.4 * concentration_matrix['OC']
        uorganic_mass = 1.4 * uncertainty_matrix['OC']
         
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (2.2 * concentration_matrix['Al'] + 2.49 * concentration_matrix['Si'] +
                               1.63 * concentration_matrix['Ca'] + 1.94 * concentration_matrix['Ti'] +
                               2.42 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                               1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                               2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals}
        
        uclosure = np.linalg.norm( [ 4.125 * uncertainty_matrix['S'],
                                    1.4 * uncertainty_matrix['OC'],
                                    uncertainty_matrix['EC'],
                                    2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                    1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                    2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        
    if equation == 'Chow_1996':
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'],
                                           uncertainty_matrix['NH4'] ], axis=0)
        
        organic_mass = 1.4 * concentration_matrix['OC']
        uorganic_mass = 1.4 * uncertainty_matrix['OC']
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (1.89 * concentration_matrix['Al'] + 2.14 * concentration_matrix['Si'] +
                               1.4 * concentration_matrix['Ca'] + 1.43 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [ 1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                                 1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'] ], axis=0)
        
        salt = concentration_matrix['Na sol'] + concentration_matrix['Na no sol']  + concentration_matrix['Cl']
        usalt = np.linalg.norm( [ uncertainty_matrix['Na sol'], uncertainty_matrix['Na no sol'], uncertainty_matrix['Cl'] ], axis=0)
        
        trace_elements = (concentration_matrix['Ti'] + concentration_matrix['V'] + concentration_matrix['Cr'] +
                          concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                          concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                          concentration_matrix['Sr'] + concentration_matrix['Pb'] + concentration_matrix['Hg'] +
                          concentration_matrix['Sb'])
        utrace_elements = np.linalg.norm( [ uncertainty_matrix['Ti'], uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                          uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                          uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                          uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                          uncertainty_matrix['Sb'] ], axis=0)
       
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'salt': salt, 'trace_elements': trace_elements}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                      'usalt': usalt, 'utrace_elements': utrace_elements}
        
        uclosure = np.linalg.norm( [uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'],
                                    1.4 * uncertainty_matrix['OC'],
                                    1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                    1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'],
                                    uncertainty_matrix['Na sol'], uncertainty_matrix['Na no sol'], uncertainty_matrix['Cl'],
                                    uncertainty_matrix['Ti'], uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                                    uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                    uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                    uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                                    uncertainty_matrix['Sb'] ], axis=0)
        
    if equation == 'Andrews_2000':
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'] ], axis=0)
        
        organic_mass = 1.4 * concentration_matrix['OC']
        uorganic_mass = 1.4 * uncertainty_matrix['OC']
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (1.89 * concentration_matrix['Al'] + 2.14 * concentration_matrix['Si'] +
                               1.4 * concentration_matrix['Ca'] + 1.43 * concentration_matrix['Fe'] +
                               1.67 * concentration_matrix['Ti'])
        ugeological_minerals = np.linalg.norm( [ 1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                                1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'],
                                                1.67 * uncertainty_matrix['Ti'] ], axis=0)
        
        trace_elements = (concentration_matrix['Cl'] + concentration_matrix['Na sol'] + concentration_matrix['Na no sol'] +
                          concentration_matrix['K'] + concentration_matrix['V'] + concentration_matrix['Cr'] +
                          concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                          concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                          concentration_matrix['Sr'] + concentration_matrix['Pb'] + concentration_matrix['Hg'] +
                          concentration_matrix['Sb'])
        utrace_elements = np.linalg.norm( [uncertainty_matrix['Cl'], uncertainty_matrix['Na sol'], uncertainty_matrix['Na no sol'], uncertainty_matrix['K'],
                                           uncertainty_matrix['V'], uncertainty_matrix['Cr'], uncertainty_matrix['Hg'], 
                                           uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                           uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                           uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Sb'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'trace_elements': trace_elements}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                      'utrace_elements': utrace_elements}
        
        uclosure = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'],
                                     1.4 * uncertainty_matrix['OC'],
                                     uncertainty_matrix['EC'],
                                     1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                     1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'],
                                     1.67 * uncertainty_matrix['Ti'],
                                     uncertainty_matrix['Cl'], uncertainty_matrix['Na sol'], uncertainty_matrix['Na no sol'], uncertainty_matrix['K'],
                                     uncertainty_matrix['V'], uncertainty_matrix['Cr'], uncertainty_matrix['Hg'], 
                                     uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                     uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                     uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Sb'] ], axis=0)
        
        
    if equation == 'Malm_2000':
        inorganic_ions = 4.125 * concentration_matrix['S'] + 1.29 * concentration_matrix['NO3']
        uinorganic_ions = np.linalg.norm( [ 4.125 * uncertainty_matrix['S'], 1.29 * uncertainty_matrix['NO3'] ], axis=0)
        
        organic_mass = 1.4 * concentration_matrix['OC']
        uorganic_mass = 1.4 * uncertainty_matrix['OC']
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (2.2 * concentration_matrix['Al'] + 2.49 * concentration_matrix['Si'] +
                               1.63 * concentration_matrix['Ca'] + 1.94 * concentration_matrix['Ti'] +
                               2.42 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [ 2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                                1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                                2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals}
    
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals}
        
        uclosure = np.linalg.norm( [ 4.125 * uncertainty_matrix['S'], 1.29 * uncertainty_matrix['NO3'],
                                    1.4 * uncertainty_matrix['OC'],
                                    uncertainty_matrix['EC'],
                                    2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                    1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                    2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        
    if equation == 'Maenhaut_2002':
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm( [uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'] ], axis=0)
        
        organic_mass = 1.4 * concentration_matrix['OC']
        uorganic_mass = 1.4 * uncertainty_matrix['OC']
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']

        geological_minerals = (2.2 * concentration_matrix['Al'] + 2.49 * concentration_matrix['Si'] +
                               1.63 * concentration_matrix['Ca'] + 1.94 * concentration_matrix['Ti'] +
                               2.42 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [ 2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                               1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                               2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        salt = concentration_matrix['Cl'] + 1.4486 * (concentration_matrix['Na sol'])
        usalt = np.linalg.norm([ uncertainty_matrix['Cl'], 1.4486 * uncertainty_matrix['Na sol']], axis=0)
        
        trace_elements = (concentration_matrix['V'] + concentration_matrix['Cr'] +
                          concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                          concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                          concentration_matrix['Sr'] + concentration_matrix['Pb'] + concentration_matrix['Hg'] +
                          concentration_matrix['Sb'])  # CHEQUEAR
        utrace_elements = np.linalg.norm( [uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                                           uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                           uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                           uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                                           uncertainty_matrix['Sb'] ], axis=0)  # CHEQUEAR
        others = concentration_matrix['K'] - 0.6 * concentration_matrix['Fe']
        uothers = np.linalg.norm( [ uncertainty_matrix['K'], 0.6 * uncertainty_matrix['Fe'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
              'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
              'salt': salt, 'trace_elements': trace_elements, 'others': others}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                      'usalt':usalt, 'utrace_elements': utrace_elements, 'uothers': uothers}
        
        uclosure = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'],
                                    1.4 * uncertainty_matrix['OC'],
                                    uncertainty_matrix['EC'],
                                    2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                    1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                    2.42 * uncertainty_matrix['Fe'],
                                    uncertainty_matrix['Cl'], 1.4486 * (uncertainty_matrix['Na sol'] + uncertainty_matrix['Na no sol']),
                                    uncertainty_matrix['K'],  uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                                    uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                    uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                    uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                                    uncertainty_matrix['Sb'] ], axis=0)
        
    if equation == 'DeBell_2006':
        inorganic_ions = 4.125 * concentration_matrix['S'] + 1.29 * concentration_matrix['NO3']
        uinorganic_ions = np.linalg.norm( [ 4.125 * uncertainty_matrix['S'], 1.29 * uncertainty_matrix['NO3'] ], axis=0)
        
        organic_mass = 1.8 * concentration_matrix['OC']
        uorganic_mass = 1.8 * uncertainty_matrix['OC']
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (2.2 * concentration_matrix['Al'] + 2.49 * concentration_matrix['Si'] +
                               1.63 * concentration_matrix['Ca'] + 1.94 * concentration_matrix['Ti'] +
                               2.42 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [ 2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                                 1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                                 2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals}
        
        uclosure = np.linalg.norm( [ 4.125 * uncertainty_matrix['S'], 1.29 * uncertainty_matrix['NO3'],
                                    1.8 * uncertainty_matrix['OC'],
                                    uncertainty_matrix['EC'],
                                    2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                    1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                    2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
            
    if equation == 'Hand_2011':
        inorganic_ions = 1.375 * concentration_matrix['SO4'] + 1.29 * concentration_matrix['NO3']
        uinorganic_ions = np.linalg.norm( [ 1.375 * uncertainty_matrix['SO4'], 1.29 * uncertainty_matrix['NO3'] ] ,axis=0)
        organic_mass = 1.8 * concentration_matrix['OC']
        uorganic_mass = 1.8 * uncertainty_matrix['OC'] # Un solo elemento
        elemental_C = concentration_matrix['EC'] # Un solo elemento
        uelemental_C = uncertainty_matrix['EC']
        geological_minerals = (2.2 * concentration_matrix['Al'] + 2.49 * concentration_matrix['Si'] +
                               1.63 * concentration_matrix['Ca'] + 1.94 * concentration_matrix['Ti'] +
                               2.42 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                               1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                               2.42 * uncertainty_matrix['Fe'] ], axis=0)
        salt = 1.8 * concentration_matrix['Cl']
        usalt = 1.8 * uncertainty_matrix['Cl'] # Un solo elemento
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'salt': salt}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C,
                       'ugeological_minerals': ugeological_minerals,
                       'usalt': usalt}
        
        uclosure = np.linalg.norm( [ 1.375 * uncertainty_matrix['SO4'], 1.29 * uncertainty_matrix['NO3'], 
                                       1.8 * uncertainty_matrix['OC'],
                                       uncertainty_matrix['EC'],
                                       2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                       1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                       2.42 * uncertainty_matrix['Fe'],
                                       1.8 * uncertainty_matrix['Cl'] ], axis=0)
        
    if equation == 'Hand_2011_mod':
        inorganic_ions = 1.375 * concentration_matrix['SO4'] + 1.29 * concentration_matrix['NO3']
        uinorganic_ions = np.linalg.norm( [ 1.375 * uncertainty_matrix['SO4'], 1.29 * uncertainty_matrix['NO3'] ] ,axis=0)
        organic_mass = 1.8 * concentration_matrix['OC']
        uorganic_mass = 1.8 * uncertainty_matrix['OC'] # Un solo elemento
        elemental_C = concentration_matrix['EC'] # Un solo elemento
        uelemental_C = uncertainty_matrix['EC']
        geological_minerals = (2.2 * concentration_matrix['Al'] + 2.49 * concentration_matrix['Si'] +
                               1.63 * concentration_matrix['Ca'] + 1.94 * concentration_matrix['Ti'] +
                               2.42 * concentration_matrix['Fe'] + 1.35 * concentration_matrix['Na no sol'])
        ugeological_minerals = np.linalg.norm( [2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                               1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                               2.42 * uncertainty_matrix['Fe'], 1.35 * uncertainty_matrix['Na no sol'] ], axis=0)
        salt = 1.8 * concentration_matrix['Cl']
        usalt = 1.8 * uncertainty_matrix['Cl'] # Un solo elemento
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'salt': salt}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C,
                       'ugeological_minerals': ugeological_minerals,
                       'usalt': usalt}
        
        uclosure = np.linalg.norm( [ 1.375 * uncertainty_matrix['SO4'], 1.29 * uncertainty_matrix['NO3'], 
                                       1.8 * uncertainty_matrix['OC'],
                                       uncertainty_matrix['EC'],
                                       2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                       1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                       2.42 * uncertainty_matrix['Fe'],
                                       1.8 * uncertainty_matrix['Cl'] , 
                                       1.35 * uncertainty_matrix['Na no sol']], axis=0)
               
    if equation == 'Simon_2011':
        # Consider all SO4 and NO3 with ammonium as counterion
        
        if "(NH4)2SO4" not in concentration_matrix:
            # Assume all SO4 is (NH4)2SO4
            concentration_matrix["(NH4)2SO4"] = (132.14 / 96.06 ) * concentration_matrix["SO4"]
            uncertainty_matrix["(NH4)2SO4"] = (132.14 / 96.06 ) * uncertainty_matrix["SO4"]
        if "NH4NO3" not in concentration_matrix:
            # Assume all NO3 is NH4NO3
            concentration_matrix["NH4NO3"] =  (80.043 / 62.004 ) * concentration_matrix["NO3"]
            uncertainty_matrix["NH4NO3"] = (80.043 / 62.004 ) * uncertainty_matrix["NO3"]
        
        inorganic_ions = concentration_matrix['(NH4)2SO4'] + concentration_matrix['NH4NO3']
        uinorganic_ions = np.linalg.norm( [ uncertainty_matrix['(NH4)2SO4'], uncertainty_matrix['NH4NO3'] ], axis=0)
        
        organic_mass = (1.8 * concentration_matrix['OC'] )
        uorganic_mass = (1.8 * uncertainty_matrix['OC'] )
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (3.48 * concentration_matrix['Si'] + 1.63 * concentration_matrix['Ca'] +
                               2.42 * concentration_matrix['Fe'] + 1.94 * concentration_matrix['Ti'])
        ugeological_minerals = np.linalg.norm( [3.48 * uncertainty_matrix['Si'], 1.63 * uncertainty_matrix['Ca'],
                               2.42 * uncertainty_matrix['Fe'], 1.94 * uncertainty_matrix['Ti'] ], axis=0)
        
        salt = 1.8 * concentration_matrix['Cl']
        usalt = 1.8 * uncertainty_matrix['Cl']
        
        others = 1.2 * (concentration_matrix['K'] - 0.6 * concentration_matrix['Fe'])
        uothers = np.linalg.norm( [1.2 * uncertainty_matrix['K'], 1.2 * 0.6 * uncertainty_matrix['Fe'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
              'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
              'salt': salt, 'others': others}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 'usalt': usalt, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals, 'uothers': uothers}
        
        uclosure = np.linalg.norm( [ uncertainty_matrix['(NH4)2SO4'], uncertainty_matrix['NH4NO3'],
                                    1.8 * uncertainty_matrix['OC'],
                                    uncertainty_matrix['EC'],
                                    3.48 * uncertainty_matrix['Si'], 1.63 * uncertainty_matrix['Ca'],
                                    2.42 * uncertainty_matrix['Fe'], 1.94 * uncertainty_matrix['Ti'],
                                    1.8 * uncertainty_matrix['Cl'],
                                    1.2 * uncertainty_matrix['K'], 1.2 * 0.6 * uncertainty_matrix['Fe'] ], axis=0)
        
        
    closure = sum(categories.values())
    # categories['unexplained'] = concentration_matrix['PM2.5'] - closure
    # ucategories['uunexplained'] = np.linalg.norm( [uclosure, uncertainty_matrix['PM2.5'] ], axis=0)
    return closure, categories, uclosure, ucategories



def mass_reconstruction_mod(conc_matrix, unc_matrix, events, equation='Simon_2011', event_labels=['SP', 'S', 'SN','SL'], 
                            event_column="Event", omoc_event=2.6, omoc_noevent=2, omoc_all=2.3, 
                            betas_event=[1,1,1,1], betas_noevent=[1,1,1,1], betas_all =[1,1,1,1],
                            all_together=False):
    
    """
    Reconstructs the mass using methodologies as described in Chow 2015. It requires a concentration
    and an uncertanty matrices.
    """

    closure = 0
    concentration_matrix = conc_matrix.copy()
    uncertainty_matrix = unc_matrix.copy()
#    data_df2 = data_df.fillna(0)
    if 'Si' not in concentration_matrix:
        concentration_matrix['Si'] = 2.174 * concentration_matrix["Al"] #2.4729 * concentration_matrix['Al']
        uncertainty_matrix['Si'] = 2.174 * concentration_matrix["Al"]#2.4729 * uncertainty_matrix['Al'] * 2 # Ese * 2 es solo para agrandar la incert Si
    concentration_matrix['S'] = concentration_matrix['SO4']/3 
    uncertainty_matrix['S'] = uncertainty_matrix['SO4']/3 
    concentration_matrix['Sr'] = 0
    uncertainty_matrix['Sr'] = 0
    concentration_matrix['Hg'] = 0
    uncertainty_matrix['Hg'] = 0
    
    if equation == 'Macias_1981':

                
        if "(NH4)2SO4" not in concentration_matrix:
            # Assume all SO4 is (NH4)2SO4
            concentration_matrix["(NH4)2SO4"] = (132.14 / 96.06 ) * concentration_matrix["SO4"]
            uncertainty_matrix["(NH4)2SO4"] = (132.14 / 96.06 ) * uncertainty_matrix["SO4"]
        if "NH4NO3" not in concentration_matrix:
            # Assume all NO3 is NH4NO3
            concentration_matrix["NH4NO3"] =  (80.043 / 62.004 ) * concentration_matrix["NO3"]
            uncertainty_matrix["NH4NO3"] = (80.043 / 62.004 ) * uncertainty_matrix["NO3"]
        
        
        inorganic_ions = concentration_matrix['(NH4)2SO4'] + concentration_matrix['NH4NO3']
        uinorganic_ions = np.linalg.norm( [uncertainty_matrix['(NH4)2SO4'], uncertainty_matrix['NH4NO3'] ], axis=0)
        
        if all_together == True:
            organic_mass = (omoc_all * concentration_matrix['OC'])
            uorganic_mass = (omoc_all * uncertainty_matrix['OC'])
        else:
            organic_mass = (omoc_noevent * concentration_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                            + omoc_event * concentration_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
            uorganic_mass = (omoc_noevent * uncertainty_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                        + omoc_event * uncertainty_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (1.89 * concentration_matrix['Al'] + 2.14 * concentration_matrix['Si'] +
                               1.4 * concentration_matrix['Ca'] + 1.2 * concentration_matrix['K'] +
                               1.43 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                                1.4 * uncertainty_matrix['Ca'], 1.2 * uncertainty_matrix['K'],
                                                1.43 * uncertainty_matrix['Fe'] ], axis=0)
        
        trace_elements = (1.25 * concentration_matrix['Cu'] + 1.24 * concentration_matrix['Zn'] +
                          1.08 * concentration_matrix['Pb'])
        utrace_elements = np.linalg.norm( [1.25 * uncertainty_matrix['Cu'], 1.24 * uncertainty_matrix['Zn'],
                                              1.08 * uncertainty_matrix['Pb'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                       'trace_elements': trace_elements}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                       'utrace_elements': utrace_elements}
        
        uclosure = np.linalg.norm( [ uncertainty_matrix['(NH4)2SO4'], uncertainty_matrix['NH4NO3'],
                                        uorganic_mass,
                                        uncertainty_matrix['EC'],
                                        1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                        1.4 * uncertainty_matrix['Ca'], 1.2 * uncertainty_matrix['K'],
                                        1.43 * uncertainty_matrix['Fe'],
                                        1.25 * uncertainty_matrix['Cu'], 1.24 * uncertainty_matrix['Zn'],
                                        1.08 * uncertainty_matrix['Pb'] ], axis=0)
                               
      
    if equation == 'Solomon_1989':
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm([ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'] ], axis=0)
        
        if all_together == True:
            organic_mass = (omoc_all * concentration_matrix['OC'])
            uorganic_mass = (omoc_all * uncertainty_matrix['OC'])
        else:
            organic_mass = (omoc_noevent * concentration_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                            + omoc_event * concentration_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
            uorganic_mass = (omoc_noevent * uncertainty_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                        + omoc_event * uncertainty_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (1.89 * concentration_matrix['Al'] + 2.14 * concentration_matrix['Si'] +
                               1.4 * concentration_matrix['Ca'] + 1.43 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [1.89 * uncertainty_matrix['Al'] + 2.14 * uncertainty_matrix['Si'] +
                               1.4 * uncertainty_matrix['Ca'] + 1.43 * uncertainty_matrix['Fe'] ], axis=0)
        
        # Pide que sean por XRF salvo Na+ y Mg++, que deberían ser por AAS. No es el caso en ninguno
        trace_elements = (concentration_matrix['Cl'] + concentration_matrix['Na sol'] + concentration_matrix['Na no sol'] +
                          concentration_matrix['K'] +
                          concentration_matrix['Ti'] + concentration_matrix['V'] + concentration_matrix['Cr'] +
                          concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                          concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                          concentration_matrix['Sr'] + concentration_matrix['Pb'] + concentration_matrix['Hg'] +
                          concentration_matrix['Sb'])
        utrace_elements = np.linalg.norm( [uncertainty_matrix['Cl'] + uncertainty_matrix['Na sol'] + uncertainty_matrix['Na no sol'] +
                                           uncertainty_matrix['K'] +
                          uncertainty_matrix['Ti'] + uncertainty_matrix['V'] + uncertainty_matrix['Cr'] +
                          uncertainty_matrix['Mn'] + uncertainty_matrix['Ni'] + uncertainty_matrix['Cu'] +
                          uncertainty_matrix['Zn'] + uncertainty_matrix['As'] + uncertainty_matrix['Se'] +
                          uncertainty_matrix['Sr'] + uncertainty_matrix['Pb'] + uncertainty_matrix['Hg'] +
                          uncertainty_matrix['Sb'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'trace_elements': trace_elements}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                      'utrace_elements': utrace_elements}
        
        uclosure = np.linalg.norm( [uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'],
                                    uorganic_mass,
                                    uncertainty_matrix['EC'],
                                    1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                    1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'],
                                    uncertainty_matrix['Cl'] + uncertainty_matrix['Na sol'] + uncertainty_matrix['Na no sol'] +
                                    uncertainty_matrix['K'] +
                                    uncertainty_matrix['Ti'] + uncertainty_matrix['V'] + uncertainty_matrix['Cr'] +
                                    uncertainty_matrix['Mn'] + uncertainty_matrix['Ni'] + uncertainty_matrix['Cu'] +
                                    uncertainty_matrix['Zn'] + uncertainty_matrix['As'] + uncertainty_matrix['Se'] +
                                    uncertainty_matrix['Sr'] + uncertainty_matrix['Pb'] + uncertainty_matrix['Hg'] +
                                    uncertainty_matrix['Sb'] ], axis=0)
                                    
        
    if equation == 'Chow_1994':
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'] ], axis=0)
        
        if all_together == True:
            organic_mass = (omoc_all * concentration_matrix['OC'])
            uorganic_mass = (omoc_all * uncertainty_matrix['OC'])
        else:
            organic_mass = (omoc_noevent * concentration_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                            + omoc_event * concentration_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
            uorganic_mass = (omoc_noevent * uncertainty_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                        + omoc_event * uncertainty_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (1.89 * concentration_matrix['Al'] + 2.14 * concentration_matrix['Si'] +
                               1.4 * concentration_matrix['Ca'] + 1.43 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'] ], axis=0)
        
        trace_elements = (concentration_matrix['Cl'] + concentration_matrix['Na sol'] + concentration_matrix['Na no sol'] +
                          concentration_matrix['K'] +
                          concentration_matrix['Ti'] + concentration_matrix['V'] + concentration_matrix['Cr'] +
                          concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                          concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                          concentration_matrix['Sr'] + concentration_matrix['Pb'] + concentration_matrix['Hg'] +
                          concentration_matrix['Sb'])
        utrace_elements = np.linalg.norm( [ uncertainty_matrix['Cl'], uncertainty_matrix['Na sol'],  uncertainty_matrix['Na no sol'],
                                            uncertainty_matrix['K'],
                                            uncertainty_matrix['Ti'], uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                                            uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                            uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                            uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                                            uncertainty_matrix['Sb'] ], axis=0)     
        
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'trace_elements': trace_elements}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                       'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                       'utrace_elements': utrace_elements}
        # Pablo, este uclosure era el único multiplicado por 1.4
        uclosure = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'],
                                    uorganic_mass,
                                    uncertainty_matrix['EC'],
                                    1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                    1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'],
                                    uncertainty_matrix['Cl'], uncertainty_matrix['Na sol'], uncertainty_matrix['Na no sol'],
                                    uncertainty_matrix['K'],
                                    uncertainty_matrix['Ti'], uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                                    uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                    uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                    uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                                    uncertainty_matrix['Sb'] ], axis=0)
                                    
        
    if equation == 'Malm_1994':
        inorganic_ions = 4.125 * concentration_matrix['S']
        uinorganic_ions = 4.125 * uncertainty_matrix['S']
        
        if all_together == True:
            organic_mass = (omoc_all * concentration_matrix['OC'])
            uorganic_mass = (omoc_all * uncertainty_matrix['OC'])
        else:
            organic_mass = (omoc_noevent * concentration_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                            + omoc_event * concentration_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
            uorganic_mass = (omoc_noevent * uncertainty_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                        + omoc_event * uncertainty_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
         
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (2.2 * concentration_matrix['Al'] + 2.49 * concentration_matrix['Si'] +
                               1.63 * concentration_matrix['Ca'] + 1.94 * concentration_matrix['Ti'] +
                               2.42 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                               1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                               2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals}
        
        uclosure = np.linalg.norm( [ 4.125 * uncertainty_matrix['S'],
                                    uorganic_mass,
                                    uncertainty_matrix['EC'],
                                    2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                    1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                    2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        
    if equation == 'Chow_1996':
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'],
                                           uncertainty_matrix['NH4'] ], axis=0)
        
        if all_together == True:
            organic_mass = (omoc_all * concentration_matrix['OC'])
            uorganic_mass = (omoc_all * uncertainty_matrix['OC'])
        else:
            organic_mass = (omoc_noevent * concentration_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                            + omoc_event * concentration_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
            uorganic_mass = (omoc_noevent * uncertainty_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                        + omoc_event * uncertainty_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (1.89 * concentration_matrix['Al'] + 2.14 * concentration_matrix['Si'] +
                               1.4 * concentration_matrix['Ca'] + 1.43 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [ 1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                                 1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'] ], axis=0)
        
        salt = concentration_matrix['Na sol'] + concentration_matrix['Na no sol']  + concentration_matrix['Cl']
        usalt = np.linalg.norm( [ uncertainty_matrix['Na sol'], uncertainty_matrix['Na no sol'], uncertainty_matrix['Cl'] ], axis=0)
        
        trace_elements = (concentration_matrix['Ti'] + concentration_matrix['V'] + concentration_matrix['Cr'] +
                          concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                          concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                          concentration_matrix['Sr'] + concentration_matrix['Pb'] + concentration_matrix['Hg'] +
                          concentration_matrix['Sb'])
        utrace_elements = np.linalg.norm( [ uncertainty_matrix['Ti'], uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                          uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                          uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                          uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                          uncertainty_matrix['Sb'] ], axis=0)
       
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'salt': salt, 'trace_elements': trace_elements}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                      'usalt': usalt, 'utrace_elements': utrace_elements}
        
        uclosure = np.linalg.norm( [uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'],
                                    uorganic_mass,
                                    1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                    1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'],
                                    uncertainty_matrix['Na sol'], uncertainty_matrix['Na no sol'], uncertainty_matrix['Cl'],
                                    uncertainty_matrix['Ti'], uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                                    uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                    uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                    uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                                    uncertainty_matrix['Sb'] ], axis=0)
        
    if equation == 'Andrews_2000':
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'] ], axis=0)
        
        if all_together == True:
            organic_mass = (omoc_all * concentration_matrix['OC'])
            uorganic_mass = (omoc_all * uncertainty_matrix['OC'])
        else:
            organic_mass = (omoc_noevent * concentration_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                            + omoc_event * concentration_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
            uorganic_mass = (omoc_noevent * uncertainty_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                        + omoc_event * uncertainty_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (1.89 * concentration_matrix['Al'] + 2.14 * concentration_matrix['Si'] +
                               2 * concentration_matrix['Ca'] + 1.43 * concentration_matrix['Fe'] +
                               1.67 * concentration_matrix['Ti'])
        ugeological_minerals = np.linalg.norm( [ 1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                                1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'],
                                                1.67 * uncertainty_matrix['Ti'] ], axis=0)
        
        trace_elements = (concentration_matrix['Cl'] + concentration_matrix['Na sol'] + concentration_matrix['Na no sol'] +
                          concentration_matrix['K'] + concentration_matrix['V'] + concentration_matrix['Cr'] +
                          concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                          concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                          concentration_matrix['Sr'] + concentration_matrix['Pb'] + concentration_matrix['Hg'] +
                          concentration_matrix['Sb'])
        utrace_elements = np.linalg.norm( [uncertainty_matrix['Cl'], uncertainty_matrix['Na sol'], uncertainty_matrix['Na no sol'], uncertainty_matrix['K'],
                                           uncertainty_matrix['V'], uncertainty_matrix['Cr'], uncertainty_matrix['Hg'], 
                                           uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                           uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                           uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Sb'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'trace_elements': trace_elements}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                      'utrace_elements': utrace_elements}
        
        uclosure = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'],
                                     uorganic_mass,
                                     uncertainty_matrix['EC'],
                                     1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                     1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'],
                                     1.67 * uncertainty_matrix['Ti'],
                                     uncertainty_matrix['Cl'], uncertainty_matrix['Na sol'], uncertainty_matrix['Na no sol'], uncertainty_matrix['K'],
                                     uncertainty_matrix['V'], uncertainty_matrix['Cr'], uncertainty_matrix['Hg'], 
                                     uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                     uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                     uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Sb'] ], axis=0)
        
        
    if equation == 'Malm_2000':
        inorganic_ions = 4.125 * concentration_matrix['S'] + 1.29 * concentration_matrix["NO3"]
        uinorganic_ions = np.linalg.norm( [ 4.125 * uncertainty_matrix['S'], 1.29 * uncertainty_matrix['NO3'] ], axis=0)
        
        if all_together == True:
            organic_mass = (omoc_all * concentration_matrix['OC'])
            uorganic_mass = (omoc_all * uncertainty_matrix['OC'])
        else:
            organic_mass = (omoc_noevent * concentration_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                            + omoc_event * concentration_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
            uorganic_mass = (omoc_noevent * uncertainty_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                        + omoc_event * uncertainty_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (2.2 * concentration_matrix['Al'] + 2.49 * concentration_matrix['Si'] +
                               1.63 * concentration_matrix['Ca'] + 1.94 * concentration_matrix['Ti'] +
                               2.42 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [ 2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                                1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                                2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals}
    
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals}
        
        uclosure = np.linalg.norm( [ 4.125 * uncertainty_matrix['S'], 1.29 * uncertainty_matrix['NO3'],
                                    uorganic_mass,
                                    uncertainty_matrix['EC'],
                                    2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                    1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                    2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        
    if (equation == 'Maenhaut_2002' or equation=='Maenhaut_2002_mod'):
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm( [uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'] ], axis=0)
        
        if all_together == True:
            organic_mass = (omoc_all * concentration_matrix['OC'])
            uorganic_mass = (omoc_all * uncertainty_matrix['OC'])
        else:
            organic_mass = (omoc_noevent * concentration_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                            + omoc_event * concentration_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
            uorganic_mass = (omoc_noevent * uncertainty_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                        + omoc_event * uncertainty_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']

        geological_minerals = (2.2 * concentration_matrix['Al'] + 2.49 * concentration_matrix['Si'] +
                               1.63 * concentration_matrix['Ca'] + 1.94 * concentration_matrix['Ti'] +
                               2.42 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [ 2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                               1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                               2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        salt = concentration_matrix['Cl'] + 1.4486 * (concentration_matrix['Na sol'])
        usalt = np.linalg.norm([ uncertainty_matrix['Cl'], 1.4486 * uncertainty_matrix['Na sol']], axis=0)
        
        trace_elements = (concentration_matrix['V'] + concentration_matrix['Cr'] +
                          concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                          concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                          concentration_matrix['Sr'] + concentration_matrix['Pb'] + concentration_matrix['Hg'] +
                          concentration_matrix['Sb'])  # CHEQUEAR
        utrace_elements = np.linalg.norm( [uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                                           uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                           uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                           uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                                           uncertainty_matrix['Sb'] ], axis=0)  # CHEQUEAR
        others = concentration_matrix['K'] - 0.6 * concentration_matrix['Fe']
        uothers = np.linalg.norm( [ uncertainty_matrix['K'], 0.6 * uncertainty_matrix['Fe'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
              'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
              'salt': salt, 'trace_elements': trace_elements, 'others': others}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                      'usalt':usalt, 'utrace_elements': utrace_elements, 'uothers': uothers}
        
        uclosure = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'],
                                    uorganic_mass,
                                    uncertainty_matrix['EC'],
                                    2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                    1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                    2.42 * uncertainty_matrix['Fe'],
                                    uncertainty_matrix['Cl'], 1.4486 * (uncertainty_matrix['Na sol']),
                                    uncertainty_matrix['K'],  uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                                    uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                    uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                    uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                                    uncertainty_matrix['Sb'] ], axis=0)
        
    if equation == 'DeBell_2006':
        inorganic_ions = 4.125 * concentration_matrix['S'] + 1.29 * concentration_matrix['NO3']
        uinorganic_ions = np.linalg.norm( [ 4.125 * uncertainty_matrix['S'], 1.29 * uncertainty_matrix['NO3'] ], axis=0)
        
        if all_together == True:
            organic_mass = (omoc_all * concentration_matrix['OC'])
            uorganic_mass = (omoc_all * uncertainty_matrix['OC'])
        else:
            organic_mass = (omoc_noevent * concentration_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                            + omoc_event * concentration_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
            uorganic_mass = (omoc_noevent * uncertainty_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                        + omoc_event * uncertainty_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (2.2 * concentration_matrix['Al'] + 2.49 * concentration_matrix['Si'] +
                               1.63 * concentration_matrix['Ca'] + 1.94 * concentration_matrix['Ti'] +
                               2.42 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [ 2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                                 1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                                 2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals}
        
        uclosure = np.linalg.norm( [ 4.125 * uncertainty_matrix['S'], 1.29 * uncertainty_matrix['NO3'],
                                    uorganic_mass,
                                    uncertainty_matrix['EC'],
                                    2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                    1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                    2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
            
    if (equation == 'Hand_2011' or equation== "Hand_2011_mod"):
        inorganic_ions = 1.375 * concentration_matrix['SO4'] + 1.29 * concentration_matrix['NO3']
        uinorganic_ions = np.linalg.norm( [ 1.375 * uncertainty_matrix['SO4'], 1.29 * uncertainty_matrix['NO3'] ] ,axis=0)

        if all_together == True:
            organic_mass = (omoc_all * concentration_matrix['OC'])
            uorganic_mass = (omoc_all * uncertainty_matrix['OC'])
        else:
            organic_mass = (omoc_noevent * concentration_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                            + omoc_event * concentration_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
            uorganic_mass = (omoc_noevent * uncertainty_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                        + omoc_event * uncertainty_matrix['OC'].where(events[event_column].isin(event_labels), other=0))

        elemental_C = concentration_matrix['EC'] # Un solo elemento
        uelemental_C = uncertainty_matrix['EC']
        geological_minerals = (2.2 * concentration_matrix['Al'] + 2.49 * concentration_matrix['Si'] +
                               1.63 * concentration_matrix['Ca'] + 1.94 * concentration_matrix['Ti'] +
                               2.42 * concentration_matrix['Fe'])
        if equation == "Hand_2011_mod":
            geological_minerals = geological_minerals + concentration_matrix['Na no sol']
        if equation =="Hand_2011":    
            ugeological_minerals = np.linalg.norm( [2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                2.42 * uncertainty_matrix['Fe'] ], axis=0)
        else:
            ugeological_minerals = np.linalg.norm( [2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                2.42 * uncertainty_matrix['Fe'] , uncertainty_matrix['Na no sol']], axis=0)
        salt = 1.8 * concentration_matrix['Cl']
        usalt = 1.8 * uncertainty_matrix['Cl'] # Un solo elemento
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
                      'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
                      'salt': salt}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C,
                       'ugeological_minerals': ugeological_minerals,
                       'usalt': usalt}
        
        if equation=="Hand_2011":
            uclosure = np.linalg.norm( [ 1.375 * uncertainty_matrix['SO4'], 1.29 * uncertainty_matrix['NO3'], 
                                       uorganic_mass,
                                       uncertainty_matrix['EC'],
                                       2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                       1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                       2.42 * uncertainty_matrix['Fe'],
                                       1.8 * uncertainty_matrix['Cl'] ], axis=0)
        else:
            uclosure = np.linalg.norm( [ 1.375 * uncertainty_matrix['SO4'], 1.29 * uncertainty_matrix['NO3'], 
                                       uorganic_mass,
                                       uncertainty_matrix['EC'],
                                       2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                       1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                       2.42 * uncertainty_matrix['Fe'],
                                       1.8 * uncertainty_matrix['Cl'], uncertainty_matrix['Na no sol'] ], 
                                      axis=0)
       
    if equation == 'Simon_2011':
        # Consider all SO4 and NO3 with ammonium as counterion
        
        if "(NH4)2SO4" not in concentration_matrix:
            # Assume all SO4 is (NH4)2SO4
            concentration_matrix["(NH4)2SO4"] = (132.14 / 96.06 ) * concentration_matrix["SO4"]
            uncertainty_matrix["(NH4)2SO4"] = (132.14 / 96.06 ) * uncertainty_matrix["SO4"]
        if "NH4NO3" not in concentration_matrix:
            # Assume all NO3 is NH4NO3
            concentration_matrix["NH4NO3"] =  (80.043 / 62.004 ) * concentration_matrix["NO3"]
            uncertainty_matrix["NH4NO3"] = (80.043 / 62.004 ) * uncertainty_matrix["NO3"]
        
        inorganic_ions = concentration_matrix['(NH4)2SO4'] + concentration_matrix['NH4NO3']
        uinorganic_ions = np.linalg.norm( [ uncertainty_matrix['(NH4)2SO4'], uncertainty_matrix['NH4NO3'] ], axis=0)
        
        if all_together == True:
            organic_mass = (omoc_all * concentration_matrix['OC'])
            uorganic_mass = (omoc_all * uncertainty_matrix['OC'])
        else:
            organic_mass = (omoc_noevent * concentration_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                            + omoc_event * concentration_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
            uorganic_mass = (omoc_noevent * uncertainty_matrix['OC'].where(~events[event_column].isin(event_labels), other=0)
                        + omoc_event * uncertainty_matrix['OC'].where(events[event_column].isin(event_labels), other=0))
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        geological_minerals = (3.48 * concentration_matrix['Si'] + 1.63 * concentration_matrix['Ca'] +
                               2.42 * concentration_matrix['Fe'] + 1.94 * concentration_matrix['Ti'])
        ugeological_minerals = np.linalg.norm( [3.48 * uncertainty_matrix['Si'], 1.63 * uncertainty_matrix['Ca'],
                               2.42 * uncertainty_matrix['Fe'], 1.94 * uncertainty_matrix['Ti'] ], axis=0)
        
        salt = 1.8 * concentration_matrix['Cl']
        usalt = 1.8 * uncertainty_matrix['Cl']
        
        others = 1.2 * (concentration_matrix['K'] - 0.6 * concentration_matrix['Fe'])
        uothers = np.linalg.norm( [1.2 * uncertainty_matrix['K'], 1.2 * 0.6 * uncertainty_matrix['Fe'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
              'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
              'salt': salt, 'others': others}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 'usalt': usalt,
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals, 'uothers': uothers}
        
        # Sumo directo uorganic_mass a uclosure porque es el unico termino de esa sumatoria en ambos casos, si agrego terminos 
        # a la sumatoria hay que agregar cada termino por separado y agregar el if de si es alltogether o no. Se está trabajando con error en norma 2 
        uclosure = np.linalg.norm( [ uncertainty_matrix['(NH4)2SO4'], uncertainty_matrix['NH4NO3'],
                                    uorganic_mass,
                                    uncertainty_matrix['EC'],
                                    3.48 * uncertainty_matrix['Si'], 1.63 * uncertainty_matrix['Ca'],
                                    2.42 * uncertainty_matrix['Fe'], 1.94 * uncertainty_matrix['Ti'],
                                    1.8 * uncertainty_matrix['Cl'],
                                    1.2 * uncertainty_matrix['K'], 1.2 * 0.6 * uncertainty_matrix['Fe'] ], axis=0)
     
       
    if equation == 'Simon_2011_mod':
        # Consider all SO4 and NO3 with ammonium as counterion
        # Agregar split dataframe a df_event df_noevent
        
        if "(NH4)2SO4" not in concentration_matrix:
            # Assume all SO4 is (NH4)2SO4
            concentration_matrix["(NH4)2SO4"] = (132.14 / 96.06 ) * concentration_matrix["SO4"]
            uncertainty_matrix["(NH4)2SO4"] = (132.14 / 96.06 ) * uncertainty_matrix["SO4"]
        if "NH4NO3" not in concentration_matrix:
            # Assume all NO3 is NH4NO3
            concentration_matrix["NH4NO3"] =  (80.043 / 62.004 ) * concentration_matrix["NO3"]
            uncertainty_matrix["NH4NO3"] = (80.043 / 62.004 ) * uncertainty_matrix["NO3"]
        
        concentration_matrix_event=concentration_matrix.where(events[event_column].isin(event_labels), other=0)
        concentration_matrix_noevent=concentration_matrix.where(~events[event_column].isin(event_labels), other=0)
        uncertainty_matrix_event=uncertainty_matrix.where(events[event_column].isin(event_labels), other=0)
        uncertainty_matrix_noevent=uncertainty_matrix.where(~events[event_column].isin(event_labels), other=0)
        
               
        if all_together == True:
            organic_mass = (betas_all[0] * concentration_matrix['OC'])
            uorganic_mass = (betas_all[0] * uncertainty_matrix['OC'])
            
            inorganic_ions = (betas_all[1] * concentration_matrix['(NH4)2SO4'] + 
                              betas_all[2] * concentration_matrix['NH4NO3'])
            uinorganic_ions = np.linalg.norm( [  betas_all[1] * uncertainty_matrix['(NH4)2SO4'],  
                                               betas_all[2] * uncertainty_matrix['NH4NO3'] ], axis=0)
            
            geological_minerals = betas_all[3] * (3.48 * concentration_matrix['Si'] + 1.63 * concentration_matrix['Ca'] +
                                   2.42 * concentration_matrix['Fe'] + 1.94 * concentration_matrix['Ti'])
            ugeological_minerals = betas_all[3] * np.linalg.norm( [3.48 * uncertainty_matrix['Si'], 1.63 * uncertainty_matrix['Ca'],
                                                    2.42 * uncertainty_matrix['Fe'], 1.94 * uncertainty_matrix['Ti'] ], 
                                                  axis=0)
        else:
            organic_mass = (betas_noevent[0] * concentration_matrix_noevent['OC']
                            + betas_event[0] * concentration_matrix_event['OC'])
            uorganic_mass = (betas_noevent[0] * uncertainty_matrix_noevent['OC']
                        + betas_event[0] * uncertainty_matrix_event['OC'])
            
            inorganic_ions = (betas_event[1] * concentration_matrix_event['(NH4)2SO4'] + 
                              betas_event[2] * concentration_matrix_event['NH4NO3'] +
                              betas_noevent[1] * concentration_matrix_noevent['(NH4)2SO4'] + 
                              betas_noevent[2] * concentration_matrix_noevent['NH4NO3'])
            
            
            uinorganic_ions = np.linalg.norm( [(betas_event[1] * uncertainty_matrix_event['(NH4)2SO4'] +
                                                betas_event[2] * uncertainty_matrix_noevent['(NH4)2SO4']),
                                               (betas_noevent[1] * uncertainty_matrix_noevent['(NH4)2SO4'] + 
                                                betas_noevent[2] * uncertainty_matrix_noevent['NH4NO3'])], axis=0)        
        
            geological_minerals = (betas_event[3] * (3.48 * concentration_matrix_event['Si'] + 
                                                    1.63 * concentration_matrix_event['Ca'] +
                                                    2.42 * concentration_matrix_event['Fe'] + 
                                                    1.94 * concentration_matrix_event['Ti']) + 
                                   betas_noevent[3] * (3.48 * concentration_matrix_event['Si'] + 
                                                    1.63 * concentration_matrix_noevent['Ca'] +
                                                    2.42 * concentration_matrix_noevent['Fe'] + 
                                                    1.94 * concentration_matrix_noevent['Ti']) )
            ugeological_minerals = np.linalg.norm( [betas_event[3] * 3.48 * uncertainty_matrix_event['Si'], 
                                                    betas_event[3] * 1.63 * uncertainty_matrix_event['Ca'],
                                                    betas_event[3] * 2.42 * uncertainty_matrix_event['Fe'], 
                                                    betas_event[3] * 1.94 * uncertainty_matrix_event['Ti'],
                                                    betas_noevent[3] * 3.48 * uncertainty_matrix_noevent['Si'], 
                                                    betas_noevent[3] * 1.63 * uncertainty_matrix_noevent['Ca'],
                                                    betas_noevent[3] * 2.42 * uncertainty_matrix_noevent['Fe'], 
                                                    betas_noevent[3] * 1.94 * uncertainty_matrix_noevent['Ti']], axis=0)
        
        elemental_C = concentration_matrix['EC']
        uelemental_C = uncertainty_matrix['EC']
        
        salt = 1.8 * concentration_matrix['Cl']
        usalt = 1.8 * uncertainty_matrix['Cl']
        
        others = 1.2 * (concentration_matrix['K'] - 0.6 * concentration_matrix['Fe'])
        uothers = np.linalg.norm( [1.2 * uncertainty_matrix['K'], 1.2 * 0.6 * uncertainty_matrix['Fe'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
              'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
              'salt': salt, 'others': others}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 'usalt': usalt,
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals, 'uothers': uothers}
        
        # Sumo directo uorganic_mass a uclosure porque es el unico termino de esa sumatoria en ambos casos, si agrego terminos 
        # a la sumatoria hay que agregar cada termino por separado y agregar el if de si es alltogether o no. Se está trabajando con error en norma 2 
        if all_together == True:
            uclosure = np.linalg.norm( [ betas_all[1] * uncertainty_matrix['(NH4)2SO4'], 
                                        betas_all[2] * uncertainty_matrix['NH4NO3'],
                                        uorganic_mass,
                                        uncertainty_matrix['EC'],
                                        betas_all[3] * 3.48 * uncertainty_matrix['Si'], 
                                        betas_all[3] * 1.63 * uncertainty_matrix['Ca'],
                                        betas_all[3] * 2.42 * uncertainty_matrix['Fe'], 
                                        betas_all[3] * 1.94 * uncertainty_matrix['Ti'],
                                        1.8 * uncertainty_matrix['Cl'],
                                        1.2 * uncertainty_matrix['K'], 1.2 * 0.6 * uncertainty_matrix['Fe'] ], axis=0)
        else:
            uclosure = np.linalg.norm( [ betas_event[1] * uncertainty_matrix_event['(NH4)2SO4'] + betas_noevent[1] * uncertainty_matrix_noevent['(NH4)2SO4'], 
                                        betas_event[2] * uncertainty_matrix_event['NH4NO3'] + betas_noevent[2] * uncertainty_matrix_noevent['NH4NO3'],
                                        uorganic_mass,
                                        uncertainty_matrix['EC'],
                                        3.48 * (betas_event[3] * uncertainty_matrix_event['Si'] + betas_noevent[3] * uncertainty_matrix_noevent['Si']), 
                                        1.63 * (betas_event[3] * uncertainty_matrix_event['Ca'] + betas_noevent[3] * uncertainty_matrix_noevent['Ca']),
                                        2.42 * (betas_event[3] * uncertainty_matrix_event['Fe'] +  betas_noevent[3] * uncertainty_matrix_noevent['Fe']),
                                        1.94 * (betas_event[3] * uncertainty_matrix_event['Ti'] + betas_noevent[3] * uncertainty_matrix_noevent['Ti']),
                                        1.8 * uncertainty_matrix_event['Cl'],
                                        1.2 * uncertainty_matrix_event['K'], 
                                        1.2 * 0.6 * uncertainty_matrix_event['Fe'] ], axis=0)
     
    closure = sum(categories.values())
    # categories['unexplained'] = concentration_matrix['PM2.5'] - closure
    # ucategories['uunexplained'] = np.linalg.norm( [uclosure, uncertainty_matrix['PM2.5'] ], axis=0)
    return closure, categories, uclosure, ucategories
        

def percentage_with_err(val, totalval, uval, utotalval):
    """
    Calculate the percentage that val represents of totalval,
    and the uncertainty associated. Return a dictionary {perc, uperc}
    """
    perc = val/totalval * 100
    
    # Note that the 100 can be taken out
    uperc = 100 * np.linalg.norm([1/totalval * uval, val/(totalval**2) * utotalval], axis=0)
    
    return {'perc': perc, 'uperc': uperc}


def linear_estimation_om_oc(conc_matrix, method='Simon_2011', ssa_as_Na=False, display_latex=False):
    from IPython.display import display, Markdown, Latex
    """
    Calculate the OM/OC ratio based on Simon et al 2011, using
    a linear regresion to adjust functions from the method.
    Default is Simon 2011:
    PM25 = mOC OC + msulf (NH4)2SO4 + mnit NH4NO3 + bsoil SOIL
           + EC + 1.8 Cl- + 1.2 (K - 0.6 Fe) + e
           
    SOIL = 3.48 Si + 1.63 Ca + 2.42 Fe + 1.94 Ti
    """
    concentration_matrix = conc_matrix.copy()
    if "C Elemental" in concentration_matrix:
        concentration_matrix["EC"] = concentration_matrix["C Elemental"]
    if "C Orgánico" in concentration_matrix:
        concentration_matrix["OC"] = concentration_matrix['C Orgánico']
    if "Si" not in concentration_matrix:
        concentration_matrix['Si'] = 2.4729 * concentration_matrix['Al']
    if "(NH4)2SO4" not in concentration_matrix:
            # Assume all SO4 is (NH4)2SO4
        concentration_matrix["(NH4)2SO4"] = (132.14 / 96.06 ) * concentration_matrix["SO4"]
    if "NH4NO3" not in concentration_matrix:
            # Assume all NO3 is NH4NO3
        concentration_matrix["NH4NO3"] =  (80.043 / 62.004 ) * concentration_matrix["NO3"]
    
    
    if method=="Macias_1981":
        concentration_matrix = concentration_matrix.dropna(subset=["(NH4)2SO4", "NH4NO3", "OC", "Al", "Si",
                                                                   "Ca", "K", "Fe", "Cu", "Zn", "Pb", "PM2.5",
                                                                   "EC"], axis=0)
        
    
        soil = (1.89 * concentration_matrix["Al"] + 2.14 * concentration_matrix["Si"] +
                1.4 * concentration_matrix["Ca"] + 1.2 * concentration_matrix["K"] +
                1.43 * concentration_matrix["Fe"])
        
        trace = (1.25 * concentration_matrix["Cu"] + 1.24 * concentration_matrix["Zn"] +
                 1.08 * concentration_matrix["Pb"])
        
        
        intercept_base = (concentration_matrix["EC"] + trace + concentration_matrix["(NH4)2SO4"] + 
                          concentration_matrix["NH4NO3"] + soil.values)
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        x = concentration_matrix["OC"].values
        mask = ~np.isnan(x) & ~np.isnan(y)
        results = linregress(x[mask], y[mask])
           
    
    if method=="Chow_1994" or method=="Solomon_1989":
        
        
        concentration_matrix = concentration_matrix.dropna(subset=["SO4", "NO3", "NH4", "Si", "Al",
                                                                   "Ca", "Ti", "Fe", "Cl", "Na sol",
                                                                   "V", "Cr", "Mn", "Ni", 
                                                                   "Cu", "Zn", "As", "Se", 
                                                                   "Pb", "Sb", "PM2.5"], axis=0)
        soil = (1.89 * concentration_matrix["Al"] + 2.14 * concentration_matrix["Si"] +
                1.4 * concentration_matrix["Ca"] +
                1.43 * concentration_matrix["Fe"])
        trace = (concentration_matrix['V'] + concentration_matrix['Cr'] +
                 concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                 concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                 concentration_matrix['Pb'] + concentration_matrix["Na sol"] + concentration_matrix["Cl"] + 
                 concentration_matrix['Sb'])
        
        intercept_base = (concentration_matrix["EC"] + trace + concentration_matrix["SO4"] 
                          + concentration_matrix["NO3"] + concentration_matrix["NH4"] + soil.values)
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        x = concentration_matrix["OC"].values
        mask = ~np.isnan(x) & ~np.isnan(y)
        results = linregress(x[mask], y[mask])
           
    
    if method=="Malm_1994":
        concentration_matrix = concentration_matrix.dropna(subset=["(NH4)2SO4", "Si", "Al",
                                                                   "Ca", "Fe", "Ti", "EC",
                                                                   "Fe", "PM2.5", "OC"], axis=0)
        soil = (2.2 * concentration_matrix["Al"] + 2.49 * concentration_matrix["Si"] +
                1.63 * concentration_matrix["Ca"] + 1.94 * concentration_matrix["Ti"] +
                2.42 * concentration_matrix["Fe"])
        intercept_base = (concentration_matrix["EC"]+concentration_matrix["(NH4)2SO4"]+soil.values)
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        x = concentration_matrix["OC"].values
        mask = ~np.isnan(x) & ~np.isnan(y)
        results = linregress(x[mask], y[mask])    
    
    if method=="Chow_1996":
        
        concentration_matrix = concentration_matrix.dropna(subset=["SO4", "NO3", "NH4", "Si", "Al",
                                                                   "Ca", "Ti", "Fe", "Cl", "Na sol",
                                                                   "V", "Cr", "Mn", "Ni", 
                                                                   "Cu", "Zn", "As", "Se", 
                                                                   "Pb", "Sb", "PM2.5"], axis=0)
        soil = (1.89 * concentration_matrix["Al"] + 2.14 * concentration_matrix["Si"] +
                1.4 * concentration_matrix["Ca"] +
                1.43 * concentration_matrix["Fe"])
        if not ssa_as_Na:
            salt = concentration_matrix["Cl"] + concentration_matrix["Na sol"]
        else:
            salt = 2.54 * concentration_matrix["Na sol"]
        trace = (concentration_matrix['V'] + concentration_matrix['Cr'] +
                 concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                 concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                 concentration_matrix['Pb'] + 
                 concentration_matrix['Sb'])
        
        intercept_base = (concentration_matrix["EC"] + salt + trace + concentration_matrix["SO4"] + 
                          concentration_matrix["NO3"] + concentration_matrix["NH4"] + soil.values )
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        x = concentration_matrix["OC"].values
        mask = ~np.isnan(x) & ~np.isnan(y)
        results = linregress(x[mask], y[mask])
           
    if method=="Andrews_2000":
        
        concentration_matrix = concentration_matrix.dropna(subset=["SO4", "NO3", "NH4", "Si", "Al",
                                                                   "Ca", "Ti", "Fe", "Cl", "Na sol",
                                                                   "V", "Cr", "Mn", "Ni",
                                                                   "Cu", "Zn", "As", "Se", 
                                                                   "Pb", "Sb", "K", "PM2.5"], axis=0)
        soil = (1.89 * concentration_matrix["Al"] + 2.14 * concentration_matrix["Si"] +
                1.4 * concentration_matrix["Ca"] + 1.67 * concentration_matrix["Ti"] +
                1.43 * concentration_matrix["Fe"])
       
        trace = (concentration_matrix['V'] + concentration_matrix['Cr'] + concentration_matrix["Na sol"] +
                 concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                 concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                 concentration_matrix['Pb'] + concentration_matrix["Cl"] +
                 concentration_matrix['Sb'])
        
        intercept_base = (concentration_matrix["EC"]  + trace + soil.values +
                         concentration_matrix["SO4"] + concentration_matrix["NO3"] +
                         concentration_matrix["NH4"])
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        x = concentration_matrix["OC"].values
        mask = ~np.isnan(x) & ~np.isnan(y)
        results = linregress(x[mask], y[mask])
         
    if (method=="Maenhaut_2002" or method=="Maenhaut_2002_mod"): # Por algún motivo tuve que sacar Sr y Hg
        concentration_matrix = concentration_matrix.dropna(subset=["SO4", "NO3", "NH4", "Si", "Al",
                                                                   "Ca", "Ti", "Fe", "Cl", "Na sol",
                                                                   "V", "Cr", "Mn", "Ni", 
                                                                   "Cu", "Zn", "As", "Se", 
                                                                   "Pb", "Sb", "K", "PM2.5"], axis=0)
        soil = (2.2 * concentration_matrix["Al"] + 2.49 * concentration_matrix["Si"] +
                1.63 * concentration_matrix["Ca"] + 1.94 * concentration_matrix["Ti"] +
                2.42 * concentration_matrix["Fe"])
        if not ssa_as_Na:
            salt = concentration_matrix["Cl"] + 1.4486 * concentration_matrix["Na sol"]
        else:
            salt = 2.54 * concentration_matrix["Na sol"]
        trace = (concentration_matrix['V'] + concentration_matrix['Cr'] +
                 concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                 concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                 concentration_matrix['Pb'] + 
                 concentration_matrix['Sb'])
        knon = concentration_matrix["K"] - 0.6 * concentration_matrix["Fe"]
        
        intercept_base = (concentration_matrix["EC"] + salt + knon + trace + soil.values +
                         concentration_matrix["SO4"] + concentration_matrix["NO3"] + 
                         concentration_matrix["NH4"])
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        x = concentration_matrix["OC"].values
  
        mask = ~np.isnan(x) & ~np.isnan(y)
        results = linregress(x[mask], y[mask])
        
    if (method=="DeBell_2006" or method =="Malm_2000"):
        ############################
        #Using the fact that the equation actually considers all S as (NH4)2SO4 
        concentration_matrix = concentration_matrix.dropna(subset=["(NH4)2SO4", "NH4NO3", "Si", "Al",
                                                                   "Ca", "Fe", "Ti", "EC",
                                                                   "Fe", "PM2.5", "OC"], axis=0)
        soil = (2.2 * concentration_matrix["Al"] + 2.49 * concentration_matrix["Si"] +
                1.63 * concentration_matrix["Ca"] + 1.94 * concentration_matrix["Ti"] +
                2.42 * concentration_matrix["Fe"])
        intercept_base = (concentration_matrix["EC"] + concentration_matrix["(NH4)2SO4"] +
                            concentration_matrix["NH4NO3"] + soil.values )
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        x = concentration_matrix["OC"].values
        mask = ~np.isnan(x) & ~np.isnan(y)
        results = linregress(x[mask], y[mask])
    
           
    if (method=="Hand_2011" or method=="Hand_2011_mod"):
        concentration_matrix = concentration_matrix.dropna(subset=["(NH4)2SO4", "NH4NO3", "Si", "Al",
                                                                   "Ca", "Fe", "Ti", "EC", "Cl",
                                                                   "Fe", "PM2.5", "OC"], axis=0)
        soil = (2.2 * concentration_matrix["Al"] + 2.49 * concentration_matrix["Si"] +
                1.63 * concentration_matrix["Ca"] + 1.94 * concentration_matrix["Ti"] +
                2.42 * concentration_matrix["Fe"])
        
        if method=="Hand_2011_mod":
            soil= soil + concentration_matrix["Na no sol"]
        
        if not ssa_as_Na:
            salt = 1.8 * concentration_matrix["Cl"]
        else:
            salt = 2.54 * concentration_matrix["Na sol"]
        intercept_base = (concentration_matrix["EC"] + salt +concentration_matrix["(NH4)2SO4"] 
                          + concentration_matrix["NH4NO3"] + soil.values)
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        x = concentration_matrix["OC"].values
        
        mask = ~np.isnan(x) & ~np.isnan(y)
        results = linregress(x[mask], y[mask])
        
    if method=="Simon_2011":
        concentration_matrix = concentration_matrix.dropna(subset=["(NH4)2SO4", "NH4NO3", "Si",
                                                                   "Ca", "Fe", "Ti", "EC", "Cl",
                                                                   "K", "Fe", "PM2.5", "OC"], axis=0)
        

        soil = (3.48 * concentration_matrix["Si"] + 1.63 * concentration_matrix["Ca"] +
                2.42 * concentration_matrix["Fe"] + 1.94 * concentration_matrix["Ti"])
        
        if not ssa_as_Na:
            salt = 1.8 * concentration_matrix["Cl"]
        else:
            salt = 2.54 * concentration_matrix["Na sol"]
        intercept_base = (concentration_matrix["EC"] + salt +
                          1.2 * (concentration_matrix["K"] - 0.6 * concentration_matrix["Fe"]))
        
    #    print(concentration_matrix)
                
        y = (concentration_matrix['PM2.5'] - intercept_base -
            concentration_matrix["(NH4)2SO4"] -
            concentration_matrix["NH4NO3"] -
            soil.values
            ).values
        x = concentration_matrix["OC"].values
        mask = ~np.isnan(x) & ~np.isnan(y)
        results = linregress(x[mask], y[mask])

    return(results)

def estimation_om_oc(conc_matrix, method='Simon_2011', ssa_as_Na=False, display_latex=False):
    from IPython.display import display, Markdown, Latex
    """
    Calculate the OM/OC ratio based on Simon et al 2011, using
    a multiple regresion with ordinary least squares to adjust
    functions from the method.
    Default is Simon 2011:
    PM25 = mOC OC + msulf (NH4)2SO4 + mnit NH4NO3 + bsoil SOIL
           + EC + 1.8 Cl- + 1.2 (K - 0.6 Fe) + e
           
    SOIL = 3.48 Si + 1.63 Ca + 2.42 Fe + 1.94 Ti
    """
    concentration_matrix = conc_matrix.copy()
    if "C Elemental" in concentration_matrix:
        concentration_matrix["EC"] = concentration_matrix["C Elemental"]
    if "C Orgánico" in concentration_matrix:
        concentration_matrix["OC"] = concentration_matrix['C Orgánico']
    if "Si" not in concentration_matrix:
        concentration_matrix['Si'] = 2.4729 * concentration_matrix['Al']
    if "(NH4)2SO4" not in concentration_matrix:
            # Assume all SO4 is (NH4)2SO4
        concentration_matrix["(NH4)2SO4"] = (132.14 / 96.06 ) * concentration_matrix["SO4"]
    if "NH4NO3" not in concentration_matrix:
            # Assume all NO3 is NH4NO3
        concentration_matrix["NH4NO3"] =  (80.043 / 62.004 ) * concentration_matrix["NO3"]
    
    
    if method=="Macias_1981":
        concentration_matrix = concentration_matrix.dropna(subset=["(NH4)2SO4", "NH4NO3", "OC", "Al", "Si",
                                                                   "Ca", "K", "Fe", "Cu", "Zn", "Pb", "PM2.5",
                                                                   "EC"], axis=0)
        
    
        soil = (1.89 * concentration_matrix["Al"] + 2.14 * concentration_matrix["Si"] +
                1.4 * concentration_matrix["Ca"] + 1.2 * concentration_matrix["K"] +
                1.43 * concentration_matrix["Fe"])
        
        trace = (1.25 * concentration_matrix["Cu"] + 1.24 * concentration_matrix["Zn"] +
                 1.08 * concentration_matrix["Pb"])
        
        
        intercept_base = concentration_matrix["EC"] + trace
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        X = np.column_stack((concentration_matrix["OC"].values,
                             concentration_matrix["(NH4)2SO4"].values,
                             concentration_matrix["NH4NO3"].values,
                             soil.values))
        X = sm.add_constant(X)
        model = sm.OLS(y, X)
        results = model.fit()       
        if display_latex:
            print(method)
            display(Latex(f'PM$_{{{2.5}}}$ = {results.params[0].round(2)} µg m$^{{{-3}}}$ +\
                    {results.params[1].round(2)} OC +\
                    {results.params[2].round(2)} (NH$_4$)$_2$SO$_4$ +\
                    {results.params[3].round(2)} NH$_4$NO$_3$ +\
                    {results.params[4].round(2)} SOIL + trace +\
                    EC'))   
           
    
    if method=="Chow_1994" or method=="Solomon_1989":
        
        
        concentration_matrix = concentration_matrix.dropna(subset=["SO4", "NO3", "NH4", "Si", "Al",
                                                                   "Ca", "Ti", "Fe", "Cl", "Na sol",
                                                                   "V", "Cr", "Mn", "Ni", 
                                                                   "Cu", "Zn", "As", "Se", 
                                                                   "Pb", "Sb", "PM2.5"], axis=0)
        soil = (1.89 * concentration_matrix["Al"] + 2.14 * concentration_matrix["Si"] +
                1.4 * concentration_matrix["Ca"] +
                1.43 * concentration_matrix["Fe"])
        trace = (concentration_matrix['V'] + concentration_matrix['Cr'] +
                 concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                 concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                 concentration_matrix['Pb'] + concentration_matrix["Na sol"] + concentration_matrix["Cl"] + 
                 concentration_matrix['Sb'])
        
        intercept_base = concentration_matrix["EC"] + trace
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        X = np.column_stack((concentration_matrix["OC"].values,
                             concentration_matrix["SO4"].values,
                             concentration_matrix["NO3"].values,
                             concentration_matrix["NH4"].values,
                             soil.values))
        X = sm.add_constant(X)
        model = sm.OLS(y, X)
        results = model.fit()
        if display_latex:
            print(method)
            display(Latex(f'PM$_{{{2.5}}}$ = {results.params[0].round(2)} µg m$^{{{-3}}}$ +\
                    {results.params[1].round(2)} OC +\
                    {results.params[2].round(2)} SO$_4$ +\
                    {results.params[3].round(2)} NO$_3$ +\
                    {results.params[4].round(2)} NH$_4$ +\
                    {results.params[5].round(2)} SOIL + trace +\
                    EC'))   
           
    
    if method=="Malm_1994":
        concentration_matrix = concentration_matrix.dropna(subset=["(NH4)2SO4", "Si", "Al",
                                                                   "Ca", "Fe", "Ti", "EC",
                                                                   "Fe", "PM2.5", "OC"], axis=0)
        soil = (2.2 * concentration_matrix["Al"] + 2.49 * concentration_matrix["Si"] +
                1.63 * concentration_matrix["Ca"] + 1.94 * concentration_matrix["Ti"] +
                2.42 * concentration_matrix["Fe"])
        intercept_base = (concentration_matrix["EC"])
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        X = np.column_stack((concentration_matrix["OC"].values,
                            concentration_matrix["(NH4)2SO4"].values,
                            soil.values))
        X = sm.add_constant(X)
    #    print(X)
    #    print(y)
        model = sm.OLS(y, X)
        results = model.fit()
        if display_latex:
            print(method)
            display(Latex(f'PM$_{{{2.5}}}$ = {results.params[0].round(2)} µg m$^{{{-3}}}$ +\
                    {results.params[1].round(2)} OC +\
                    {results.params[2].round(2)} (NH$_4$)$_2$SO$_4$ +\
                    {results.params[3].round(2)} SOIL'))   
    
    
    if method=="Chow_1996":
        
        concentration_matrix = concentration_matrix.dropna(subset=["SO4", "NO3", "NH4", "Si", "Al",
                                                                   "Ca", "Ti", "Fe", "Cl", "Na sol",
                                                                   "V", "Cr", "Mn", "Ni", 
                                                                   "Cu", "Zn", "As", "Se", 
                                                                   "Pb", "Sb", "PM2.5"], axis=0)
        soil = (1.89 * concentration_matrix["Al"] + 2.14 * concentration_matrix["Si"] +
                1.4 * concentration_matrix["Ca"] +
                1.43 * concentration_matrix["Fe"])
        if not ssa_as_Na:
            salt = concentration_matrix["Cl"] + concentration_matrix["Na sol"]
        else:
            salt = 2.54 * concentration_matrix["Na sol"]
        trace = (concentration_matrix['V'] + concentration_matrix['Cr'] +
                 concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                 concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                 concentration_matrix['Pb'] + 
                 concentration_matrix['Sb'])
        
        intercept_base = concentration_matrix["EC"] + salt + trace
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        X = np.column_stack((concentration_matrix["OC"].values,
                             concentration_matrix["SO4"].values,
                             concentration_matrix["NO3"].values,
                             concentration_matrix["NH4"].values,
                             soil.values))
        X = sm.add_constant(X)
        model = sm.OLS(y, X)
        results = model.fit()
        if display_latex:
            print(method)
            display(Latex(f'PM$_{{{2.5}}}$ = {results.params[0].round(2)} µg m$^{{{-3}}}$ +\
                    {results.params[1].round(2)} OC +\
                    {results.params[2].round(2)} SO$_4$ +\
                    {results.params[3].round(2)} NO$_3$ +\
                    {results.params[4].round(2)} NH$_4$ +\
                    {results.params[5].round(2)} SOIL + trace +\
                    EC + Cl + Na'))
           
    
    if method=="Andrews_2000":
        
        concentration_matrix = concentration_matrix.dropna(subset=["SO4", "NO3", "NH4", "Si", "Al",
                                                                   "Ca", "Ti", "Fe", "Cl", "Na sol",
                                                                   "V", "Cr", "Mn", "Ni",
                                                                   "Cu", "Zn", "As", "Se", 
                                                                   "Pb", "Sb", "K", "PM2.5"], axis=0)
        soil = (1.89 * concentration_matrix["Al"] + 2.14 * concentration_matrix["Si"] +
                1.4 * concentration_matrix["Ca"] + 1.67 * concentration_matrix["Ti"] +
                1.43 * concentration_matrix["Fe"])
       
        trace = (concentration_matrix['V'] + concentration_matrix['Cr'] + concentration_matrix["Na sol"] +
                 concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                 concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                 concentration_matrix['Pb'] + concentration_matrix["Cl"] +
                 concentration_matrix['Sb'])
        
        intercept_base = concentration_matrix["EC"]  + trace
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        X = np.column_stack((concentration_matrix["OC"].values,
                             concentration_matrix["SO4"].values,
                             concentration_matrix["NO3"].values,
                             concentration_matrix["NH4"].values,
                             soil.values))
        X = sm.add_constant(X)
        model = sm.OLS(y, X)
        results = model.fit()
        if display_latex:
            print(method)
            display(Latex(f'PM$_{{{2.5}}}$ = {results.params[0].round(2)} µg m$^{{{-3}}}$ +\
                    {results.params[1].round(2)} OC +\
                    {results.params[2].round(2)} SO$_4$ +\
                    {results.params[3].round(2)} NO$_3$ +\
                    {results.params[4].round(2)} NH$_4$ +\
                    {results.params[5].round(2)} SOIL + trace +\
                    EC'))   
        
        
    
    if (method=="Maenhaut_2002" or method=="Maenhaut_2002_mod"): # Por algún motivo tuve que sacar Sr y Hg
        concentration_matrix = concentration_matrix.dropna(subset=["SO4", "NO3", "NH4", "Si", "Al",
                                                                   "Ca", "Ti", "Fe", "Cl", "Na sol",
                                                                   "V", "Cr", "Mn", "Ni", 
                                                                   "Cu", "Zn", "As", "Se", 
                                                                   "Pb", "Sb", "K", "PM2.5"], axis=0)
        soil = (2.2 * concentration_matrix["Al"] + 2.49 * concentration_matrix["Si"] +
                1.63 * concentration_matrix["Ca"] + 1.94 * concentration_matrix["Ti"] +
                2.42 * concentration_matrix["Fe"])
        if not ssa_as_Na:
            salt = concentration_matrix["Cl"] + 1.4486 * concentration_matrix["Na sol"]
        else:
            salt = 2.54 * concentration_matrix["Na sol"]
        trace = (concentration_matrix['V'] + concentration_matrix['Cr'] +
                 concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                 concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                 concentration_matrix['Pb'] + 
                 concentration_matrix['Sb'])
        knon = concentration_matrix["K"] - 0.6 * concentration_matrix["Fe"]
        
        intercept_base = concentration_matrix["EC"] + salt + knon + trace
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        X = np.column_stack((concentration_matrix["OC"].values,
                             concentration_matrix["SO4"].values,
                             concentration_matrix["NO3"].values,
                             concentration_matrix["NH4"].values,
                             soil.values))
        X = sm.add_constant(X)
        model = sm.OLS(y, X)
        results = model.fit()       
        if display_latex:
            print(method)
            display(Latex(f'PM$_{{{2.5}}}$ = {results.params[0].round(2)} µg m$^{{{-3}}}$ +\
                    {results.params[1].round(2)} OC +\
                    {results.params[2].round(2)} SO$_4$ +\
                    {results.params[3].round(2)} NO$_3$ +\
                    {results.params[4].round(2)} NH$_4$ +\
                    {results.params[5].round(2)} SOIL + trace + KNON +\
                    EC + Cl + 1.4486 Na'))   
        
        
    if (method=="DeBell_2006" or method =="Malm_2000"):
        ############################
        #Using the fact that the equation actually considers all S as (NH4)2SO4 
        concentration_matrix = concentration_matrix.dropna(subset=["(NH4)2SO4", "NH4NO3", "Si", "Al",
                                                                   "Ca", "Fe", "Ti", "EC",
                                                                   "Fe", "PM2.5", "OC"], axis=0)
        soil = (2.2 * concentration_matrix["Al"] + 2.49 * concentration_matrix["Si"] +
                1.63 * concentration_matrix["Ca"] + 1.94 * concentration_matrix["Ti"] +
                2.42 * concentration_matrix["Fe"])
        intercept_base = (concentration_matrix["EC"])
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        X = np.column_stack((concentration_matrix["OC"].values,
                            concentration_matrix["(NH4)2SO4"].values,
                            concentration_matrix["NH4NO3"].values,
                            soil.values))
        X = sm.add_constant(X)
    #    print(X)
    #    print(y)
        model = sm.OLS(y, X)
        results = model.fit()

        if display_latex:
            print(method)
            display(Latex(f'PM$_{{{2.5}}}$ = {results.params[0].round(2)} µg m$^{{{-3}}}$ +\
                    {results.params[1].round(2)} OC +\
                    {results.params[2].round(2)} (NH$_4$)$_2$SO$_4$ +\
                    {results.params[3].round(2)} NH$_4$NO$_3$ +\
                    {results.params[4].round(2)} SOIL'))       
           
    if (method=="Hand_2011" or method=="Hand_2011_mod"):
        concentration_matrix = concentration_matrix.dropna(subset=["(NH4)2SO4", "NH4NO3", "Si", "Al",
                                                                   "Ca", "Fe", "Ti", "EC", "Cl",
                                                                   "Fe", "PM2.5", "OC"], axis=0)
        soil = (2.2 * concentration_matrix["Al"] + 2.49 * concentration_matrix["Si"] +
                1.63 * concentration_matrix["Ca"] + 1.94 * concentration_matrix["Ti"] +
                2.42 * concentration_matrix["Fe"])
        
        if method=="Hand_2011_mod":
            soil= soil + concentration_matrix["Na no sol"]
        
        if not ssa_as_Na:
            salt = 1.8 * concentration_matrix["Cl"]
        else:
            salt = 2.54 * concentration_matrix["Na sol"]
        intercept_base = (concentration_matrix["EC"] + salt)
        
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        X = np.column_stack((concentration_matrix["OC"].values,
                            concentration_matrix["(NH4)2SO4"].values,
                            concentration_matrix["NH4NO3"].values,
                            soil.values))
        X = sm.add_constant(X)
    #    print(X)
    #    print(y)
        model = sm.OLS(y, X)
        results = model.fit()
        if display_latex:
            print(method)
            display(Latex(f'PM$_{{{2.5}}}$ = {results.params[0].round(2)} µg m$^{{{-3}}}$ +\
                    {results.params[1].round(2)} OC +\
                    {results.params[2].round(2)} (NH$_4$)$_2$SO$_4$ +\
                    {results.params[3].round(2)} NH$_4$NO$_3$ +\
                    {results.params[4].round(2)} SOIL +\
                    1.8 Cl$^-$ + EC'))       
    
    if method=="Simon_2011":
        concentration_matrix = concentration_matrix.dropna(subset=["(NH4)2SO4", "NH4NO3", "Si",
                                                                   "Ca", "Fe", "Ti", "EC", "Cl",
                                                                   "K", "Fe", "PM2.5", "OC"], axis=0)
        

        soil = (3.48 * concentration_matrix["Si"] + 1.63 * concentration_matrix["Ca"] +
                2.42 * concentration_matrix["Fe"] + 1.94 * concentration_matrix["Ti"])
        
        if not ssa_as_Na:
            salt = 1.8 * concentration_matrix["Cl"]
        else:
            salt = 2.54 * concentration_matrix["Na sol"]
        intercept_base = (concentration_matrix["EC"] + salt +
                          1.2 * (concentration_matrix["K"] - 0.6 * concentration_matrix["Fe"]))
        
    #    print(concentration_matrix)
                
        y = (concentration_matrix['PM2.5'] - intercept_base).values
        X = np.column_stack((concentration_matrix["OC"].values,
                            concentration_matrix["(NH4)2SO4"].values,
                            concentration_matrix["NH4NO3"].values,
                            soil.values))
        X = sm.add_constant(X)
    #    print(X)
    #    print(y)
        model = sm.OLS(y, X)
        results = model.fit()

        if display_latex:
            print("method")
            display(Latex(f'PM$_{{{2.5}}}$ = {results.params[0].round(2)} µg m$^{{{-3}}}$ +\
                    {results.params[1].round(2)} OC +\
                    {results.params[2].round(2)} (NH$_4$)$_2$SO$_4$ +\
                    {results.params[3].round(2)} NH$_4$NO$_3$ +\
                    {results.params[4].round(2)} SOIL +\
                    1.8 Cl$^-$ + 1.2 KNON + EC'))
    return(results)

def calculate_seasonal(conc_matrix):
    """Calculates seasonal means, an all around mean and standard
    deviations. Adds as well the number of samples per period as 
    last row"""

    matrix = conc_matrix.copy()

    matrix['month'] = pd.DatetimeIndex(matrix.index)
    matrix['month'] = matrix['month'].dt.to_period('M')

    month_to_season = {1:'DJF', 2:'DJF', 3:'MAM', 4:'MAM', 5:'MAM', 6:'JJA',
                       7:'JJA', 8:'JJA', 9:'SON', 10:'SON', 11:'SON', 12:'DJF'}

    matrix['season'] = matrix.index.month.map(month_to_season)
    #print(matrix['season'])

    matrix_seasonal = matrix.groupby(matrix['season']).mean(numeric_only=True).transpose()



    matrix_seasonal_std = matrix.groupby(matrix['season']).std(numeric_only=True).transpose()

    # Count number of samples per season
    count = matrix.groupby(matrix['season']).agg('count')['Ag'].transpose()
    count["All_average"] = len(matrix["PM2.5"])
    count = count.rename('N')

    for key in matrix_seasonal_std.keys():
        matrix_seasonal_std = matrix_seasonal_std.rename(columns={key:f'{key}_std'})

    matrix_seasonal = matrix_seasonal.join(matrix_seasonal_std)

    matrix_seasonal['All_average'] = matrix.mean(numeric_only=True, axis=0).transpose()
    matrix_seasonal['All_std'] = matrix.std(numeric_only=True, axis=0).transpose()

    matrix_seasonal = matrix_seasonal.reindex(sorted(matrix_seasonal.columns), axis=1)
    #matrix_seasonal = matrix_seasonal.append(count)
    matrix_seasonal = pd.concat([matrix_seasonal, count.to_frame().T], axis=0)
    matrix_seasonal = matrix_seasonal[[
        "All_average", "All_std", "DJF", "DJF_std", "MAM", "MAM_std",
        "JJA", "JJA_std", "SON", "SON_std",
    ]]

    return(matrix_seasonal)


def axvlines(ax=None, xs=[0, 1], ymin=0, ymax=1, **kwargs):
    ax = ax or plt.gca()
    for x in xs:
        ax.axvline(x, ymin=ymin, ymax=ymax, **kwargs)

def select_events(df, events=events):
    return df.where(events[event_columnname].isin(event_labels))


def select_no_events(df, events=events):
    return df.where(~events[event_columnname].isin(event_labels))


def average_mass_reconstruction(mass_Hand,mass_Maenhaut,mass_Simon):
    # TypeError: loop of ufunc does not support argument 0 of type dict_values which has no callable conjugate method
    categories = {}

    for key in mass_Hand[1].keys():
        categories[key] = (mass_Simon[1][key] + mass_Hand[1][key] + mass_Maenhaut[1][key])/3
    # Hand no tiene Others y Maenhaut tiene other + trace_elements
    categories['others'] = (mass_Simon[1]['others'] + (mass_Maenhaut[1]['others'] + mass_Maenhaut[1]['trace_elements']))/3
    
    ucategories = {}
    
    for key in mass_Hand[3].keys():
        ucategories[key] = np.linalg.norm([mass_Hand[3][key], mass_Maenhaut[3][key], mass_Simon[3][key]], ord=1, axis=0)/3
    
    ucategories["uothers"] = np.linalg.norm(
        [mass_Simon[3]["uothers"], mass_Maenhaut[3]["uothers"], mass_Maenhaut[3]["utrace_elements"]],
        ord=1, axis=0)/3

    closure = sum(categories.values())
    # print(pd.DataFrame.from_dict(ucategories))
    uclosure = np.linalg.norm(pd.DataFrame.from_dict(ucategories), axis=1)
    # print(uclosure)
    return closure, categories, uclosure, ucategories

def print_stats(model, X, y):
    
    params = np.append(model.intercept_,model.coef_)
    predictions = model.predict(X)


    newX = pd.DataFrame({"Constant":np.ones(len(X))}).join(pd.DataFrame(X))
    MSE = (sum((y-predictions)**2))/(len(newX)-len(newX.columns))

    var_b = MSE*(np.linalg.inv(np.dot(newX.T,newX)).diagonal())
    sd_b = np.sqrt(var_b)
    ts_b = params/ sd_b

    p_values =[2*(1-stats.t.cdf(np.abs(i),(len(newX)-len(X[0])))) for i in ts_b]

    sd_b = np.round(sd_b,3)
    ts_b = np.round(ts_b,3)
    p_values = np.round(p_values,3)
    params = np.round(params,4)

    myDF3 = pd.DataFrame()
    myDF3["Coefficients"],myDF3["Standard Errors"],myDF3["t values"],myDF3["p-value"] = [params,sd_b,ts_b,p_values]
    print(myDF3)

class Equation:
    """Class that represent a equation for mass reconstruction"""
    def __init__(self, name, filenameParameters):
        
        self.name = name
        self.ocEvents = ocEvents
        self.ocNevents = ocNoevents