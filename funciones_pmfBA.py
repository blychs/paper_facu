#!/usr/bin/env python
# coding: utf-8

# In[ ]:
 

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
from pandas.api.types import is_numeric_dtype
import xarray as xr



def corr_matrix(df, minimum=-1, maximum=1, method='pearson'):
    '''Calcula matriz de correlación y devuelve los valores entre las corr
    máximas y mínimas'''
    correlation = df.corr(method)
    return correlation.where(correlation >= minimum).where(correlation <= maximum)


def mass_reconstruction(concentration_matrix, uncertainty_matrix, equation='Hand_2011'):
    
    """
    Reconstructs the mass using methodologies as described in Chow 2015. It requires a concentration
    and an uncertanty matrices.
    """

    closure = 0
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
        inorganic_ions = concentration_matrix['(NH4)2SO4'] + concentration_matrix['NH4NO3']
        uinorganic_ions = np.linalg.norm( [uncertainty_matrix['(NH4)2SO4'], uncertainty_matrix['NH4NO3'] ], axis=0)
        
        organic_mass = 1.5 * concentration_matrix['C Orgánico']
        uorganic_mass = 1.5 * uncertainty_matrix['C Orgánico']
        
        elemental_C = concentration_matrix['C Elemental']
        uelemental_C = uncertainty_matrix['C Elemental']
        
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
        
        uclosure = np.linalg.norm( [ uncertainty_matrix['(NH4)2SO4'], uncertainty_matrix['NH4SO3'],
                                        1.5 * uncertainty_matrix['C Orgánico'],
                                        uncertainty_matrix['C Elemental'],
                                        1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                        1.4 * uncertainty_matrix['Ca'], 1.2 * uncertainty_matrix['K'],
                                        1.43 * uncertainty_matrix['Fe'],
                                        1.25 * uncertainty_matrix['Cu'], 1.24 * uncertainty_matrix['Zn'],
                                        1.08 * uncertainty_matrix['Pb'] ], axis=0)
                               
      
    if equation == 'Solomon_1989':
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm([ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'] ], axis=0)
        
        organic_mass = 1.4 * concentration_matrix['C Orgánico']
        uorganic_mass = 1.4 * uncertainty_matrix['C Orgánico']
        
        elemental_C = concentration_matrix['C Elemental']
        uelemental_C = uncertainty_matrix['C Elemental']
        
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
                                    1.4 * uncertainty_matrix['C Orgánico'],
                                    uncertainty_matrix['C Elemental'],
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
        
        organic_mass = 1.4 * concentration_matrix['C Orgánico']
        uorganic_mass = 1.4 * uncertainty_matrix['C Orgánico']
        
        elemental_C = concentration_matrix['C Elemental']
        uelemental_C = uncertainty_matrix['C Elemental']
        
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
                                    1.4 * uncertainty_matrix['C Orgánico'],
                                    uncertainty_matrix['C Elemental'],
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
        
        organic_mass = 1.4 * concentration_matrix['C Orgánico']
        uorganic_mass = 1.4 * uncertainty_matrix['C Orgánico']
         
        elemental_C = concentration_matrix['C Elemental']
        uelemental_C = uncertainty_matrix['C Elemental']
        
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
                                    1.4 * uncertainty_matrix['C Orgánico'],
                                    uncertainty_matrix['C Elemental'],
                                    2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                    1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                    2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        
    if equation == 'Chow_1996':
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'],
                                           uncertainty_matrix['NH4'] ], axis=0)
        
        organic_mass = 1.4 * concentration_matrix['C Orgánico']
        uorganic_mass = 1.4 * uncertainty_matrix['C Orgánico']
        
        elemental_C = concentration_matrix['C Elemental']
        uelemental_C = uncertainty_matrix['C Elemental']
        
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
                                    1.4 * uncertainty_matrix['C Orgánico'],
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
        
        organic_mass = 1.4 * concentration_matrix['C Orgánico']
        uorganic_mass = 1.4 * uncertainty_matrix['C Orgánico']
        
        elemental_C = concentration_matrix['C Elemental']
        uelemental_C = uncertainty_matrix['C Elemental']
        
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
                                     1.4 * uncertainty_matrix['C Orgánico'],
                                     uncertainty_matrix['C Elemental'],
                                     1.89 * uncertainty_matrix['Al'], 2.14 * uncertainty_matrix['Si'],
                                     1.4 * uncertainty_matrix['Ca'], 1.43 * uncertainty_matrix['Fe'],
                                     1.67 * uncertainty_matrix['Ti'],
                                     uncertainty_matrix['Cl'], uncertainty_matrix['Na sol'], uncertainty_matrix['Na no sol'], uncertainty_matrix['K'],
                                     uncertainty_matrix['V'], uncertainty_matrix['Cr'], uncertainty_matrix['Hg'], 
                                     uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                     uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                     uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Sb'] ], axis=0)
        
        
    if equation == 'Malm_2000':
        inorganic_ions = 1.125 * concentration_matrix['S'] + 1.29 * concentration_matrix['NO3']
        uinorganic_ions = np.linalg.norm( [ 1.125 * uncertainty_matrix['S'], 1.29 * uncertainty_matrix['NO3'] ], axis=0)
        
        organic_mass = 1.4 * concentration_matrix['C Orgánico']
        uorganic_mass = 1.4 * uncertainty_matrix['C Orgánico']
        
        elemental_C = concentration_matrix['C Elemental']
        uelemental_C = uncertainty_matrix['C Elemental']
        
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
        
        uclosure = np.linalg.norm( [ 1.125 * uncertainty_matrix['S'], 1.29 * uncertainty_matrix['NO3'],
                                    1.4 * uncertainty_matrix['C Orgánico'],
                                    uncertainty_matrix['C Elemental'],
                                    2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                    1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                    2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        
    if equation == 'Maenhaut_2002':
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm( [uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'] ], axis=0)
        
        organic_mass = 1.4 * concentration_matrix['C Orgánico']
        uorganic_mass = 1.4 * uncertainty_matrix['C Orgánico']
        
        elemental_C = concentration_matrix['C Elemental']
        uelemental_C = uncertainty_matrix['C Elemental']

        geological_minerals = (2.2 * concentration_matrix['Al'] + 2.49 * concentration_matrix['Si'] +
                               1.63 * concentration_matrix['Ca'] + 1.94 * concentration_matrix['Ti'] +
                               2.42 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [ 2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                               1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                               2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        salt = concentration_matrix['Cl'] + 1.4486 * (concentration_matrix['Na sol'])
        usalt = np.linalg.norm([ uncertainty_matrix['Cl'], 1.4486 * uncertainty_matrix['Na sol']], axis=0)
        
        trace_elements = (concentration_matrix['K'] +  concentration_matrix['V'] + concentration_matrix['Cr'] +
                          concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                          concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                          concentration_matrix['Sr'] + concentration_matrix['Pb'] + concentration_matrix['Hg'] +
                          concentration_matrix['Sb'])  # CHEQUEAR
        utrace_elements = np.linalg.norm( [ uncertainty_matrix['K'],  uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                                           uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                           uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                           uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                                           uncertainty_matrix['Sb'] ], axis=0)  # CHEQUEAR
        others = concentration_matrix['K'] - 0.6 * concentration_matrix['Fe']
        uothers = np.linalg.norm( [ uncertainty_matrix['K'], 0.6 * uncertainty_matrix['Fe'] ])
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
              'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
              'salt': salt, 'trace_elements': trace_elements, 'others': others}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                      'usalt':usalt, 'utrace_elements': utrace_elements, 'uothers': uothers}
        
        uclosure = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'],
                                    1.4 * uncertainty_matrix['C Orgánico'],
                                    uncertainty_matrix['C Elemental'],
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
        
        organic_mass = 1.8 * concentration_matrix['C Orgánico']
        uorganic_mass = 1.8 * uncertainty_matrix['C Orgánico']
        
        elemental_C = concentration_matrix['C Elemental']
        uelemental_C = uncertainty_matrix['C Elemental']
        
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
                                    1.8 * uncertainty_matrix['C Orgánico'],
                                    uncertainty_matrix['C Elemental'],
                                    2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                    1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                    2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
            
    if equation == 'Hand_2011':
        inorganic_ions = 1.375 * concentration_matrix['SO4'] + 1.29 * concentration_matrix['NO3']
        uinorganic_ions = np.linalg.norm( [ 1.375 * uncertainty_matrix['SO4'], 1.29 * uncertainty_matrix['NO3'] ] ,axis=0)
        organic_mass = 1.8 * concentration_matrix['C Orgánico']
        uorganic_mass = 1.8 * uncertainty_matrix['C Orgánico'] # Un solo elemento
        elemental_C = concentration_matrix['C Elemental'] # Un solo elemento
        uelemental_C = uncertainty_matrix['C Elemental']
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
                                       1.8 * uncertainty_matrix['C Orgánico'],
                                       uncertainty_matrix['C Elemental'],
                                       2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                       1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                       2.42 * uncertainty_matrix['Fe'],
                                       1.8 * uncertainty_matrix['Cl'] ], axis=0)
        
       
    if equation == 'Simon_2011':
        inorganic_ions = concentration_matrix['(NH4)2SO4'] + concentration_matrix['NH4NO3']
        uinorganic_ions = np.linalg.norm( [ uncertainty_matrix['(NH4)2SO4'], uncertainty_matrix['NH4NO3'] ], axis=0)
        
        organic_mass = (1.8 * concentration_matrix['C Orgánico'] +
                        0.2 * concentration_matrix['C Orgánico'].where(events['Event']=='SN', other=0) +
                        0.2 * concentration_matrix['C Orgánico'].where(events['Event']=='SP', other=0) +
                        0.4 * concentration_matrix['C Orgánico'].where(events['Event']=='S', other=0))
        uorganic_mass = (1.8 * uncertainty_matrix['C Orgánico'] +
                         0.2 * uncertainty_matrix['C Orgánico'].where(events['Event']=='SN', other=0) +
                         0.2 * uncertainty_matrix['C Orgánico'].where(events['Event']=='SP', other=0) +
                         0.4 * uncertainty_matrix['C Orgánico'].where(events['Event']=='S', other=0)) 
        
        elemental_C = concentration_matrix['C Elemental']
        uelemental_C = uncertainty_matrix['C Elemental']
        
        geological_minerals = (3.48 * concentration_matrix['Si'] + 1.63 * concentration_matrix['Ca'] +
                               2.42 * concentration_matrix['Fe'] + 1.94 * concentration_matrix['Ti'])
        geological_minerals = np.linalg.norm( [3.48 * uncertainty_matrix['Si'], 1.63 * uncertainty_matrix['Ca'],
                               2.42 * uncertainty_matrix['Fe'], 1.94 * uncertainty_matrix['Ti'] ], axis=0)
        
        salt = 1.8 * concentration_matrix['Cl']
        salt = 1.8 * uncertainty_matrix['Cl']
        
        others = 1.2 * (concentration_matrix['K'] - 0.6 * concentration_matrix['Fe'])
        others = np.linalg.norm( [1.2 * uncertainty_matrix['K'], 1.2 * 0.6 * uncertainty_matrix['Fe'] ], axis=0)
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
              'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
              'salt': salt, 'others': others}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                      'utrace_elements': utrace_elements}
        
        uclosure = np.linalg.norm( [ uncertainty_matrix['(NH4)2SO4'], uncertainty_matrix['NH4NO3'],
                                    1.8 * uncertainty_matrix['C Orgánico'],
                                    uncertainty_matrix['C Elemental'],
                                    3.48 * uncertainty_matrix['Si'], 1.63 * uncertainty_matrix['Ca'],
                                    2.42 * uncertainty_matrix['Fe'], 1.94 * uncertainty_matrix['Ti'],
                                    1.8 * uncertainty_matrix['Cl'],
                                    1.2 * uncertainty_matrix['K'], 1.2 * 0.6 * uncertainty_matrix['Fe'] ], axis=0)
        
        
    closure = sum(categories.values())
    categories['unexplained'] = concentration_matrix['PM2.5'] - closure
    ucategories['uunexplained'] = np.linalg.norm( [uclosure, uncertainty_matrix['PM2.5'] ], axis=0)
    return closure, categories, uclosure, ucategories



def mass_reconstruction_mod(concentration_matrix, uncertainty_matrix, events, equation='Hand_2011'):
    
    """
    Reconstructs the mass using methodologies as described in Chow 2015. It requires a concentration
    and an uncertanty matrices. It adds special cases if there are events.
    """

    closure = 0
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
    
        
    if equation == 'Maenhaut_2002':
        inorganic_ions = concentration_matrix['SO4'] + concentration_matrix['NO3'] + concentration_matrix['NH4']
        uinorganic_ions = np.linalg.norm( [uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'] ], axis=0)
        
        organic_mass = 1.4 * concentration_matrix['C Orgánico']
        uorganic_mass = 1.4 * uncertainty_matrix['C Orgánico']
        
        elemental_C = concentration_matrix['C Elemental']
        uelemental_C = uncertainty_matrix['C Elemental']

        geological_minerals = (2.2 * concentration_matrix['Al'] + 2.49 * concentration_matrix['Si'] +
                               1.63 * concentration_matrix['Ca'] + 1.94 * concentration_matrix['Ti'] +
                               2.42 * concentration_matrix['Fe'])
        ugeological_minerals = np.linalg.norm( [ 2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                               1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                               2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
        salt = concentration_matrix['Cl'] + 1.4486 * (concentration_matrix['Na sol'] + concentration_matrix['Na sol'])
        usalt = np.linalg.norm([ uncertainty_matrix['Cl'], 1.4486 * uncertainty_matrix['Na sol'], 1.4486 * uncertainty_matrix['Na sol']], axis=0)
        
        trace_elements = (concentration_matrix['K'] +  concentration_matrix['V'] + concentration_matrix['Cr'] +
                          concentration_matrix['Mn'] + concentration_matrix['Ni'] + concentration_matrix['Cu'] +
                          concentration_matrix['Zn'] + concentration_matrix['As'] + concentration_matrix['Se'] +
                          concentration_matrix['Sr'] + concentration_matrix['Pb'] + concentration_matrix['Hg'] +
                          concentration_matrix['Sb'])  # CHEQUEAR
        utrace_elements = np.linalg.norm( [ uncertainty_matrix['K'],  uncertainty_matrix['V'], uncertainty_matrix['Cr'],
                                           uncertainty_matrix['Mn'], uncertainty_matrix['Ni'], uncertainty_matrix['Cu'],
                                           uncertainty_matrix['Zn'], uncertainty_matrix['As'], uncertainty_matrix['Se'],
                                           uncertainty_matrix['Sr'], uncertainty_matrix['Pb'], uncertainty_matrix['Hg'],
                                           uncertainty_matrix['Sb'] ], axis=0)  # CHEQUEAR
        others = concentration_matrix['K'] - 0.6 * concentration_matrix['Fe']
        uothers = np.linalg.norm( [ uncertainty_matrix['K'], 0.6 * uncertainty_matrix['Fe'] ])
        
        categories = {'inorganic_ions': inorganic_ions, 'organic_mass': organic_mass, 
              'elemental_C': elemental_C, 'geological_minerals': geological_minerals,
              'salt': salt, 'trace_elements': trace_elements, 'others': others}
        
        ucategories = {'uinorganic_ions': uinorganic_ions, 'uorganic_mass': uorganic_mass, 
                      'uelemental_C': uelemental_C, 'ugeological_minerals': ugeological_minerals,
                      'usalt':usalt, 'utrace_elements': utrace_elements, 'uothers': uothers}
        
        uclosure = np.linalg.norm( [ uncertainty_matrix['SO4'], uncertainty_matrix['NO3'], uncertainty_matrix['NH4'],
                                    1.4 * uncertainty_matrix['C Orgánico'],
                                    uncertainty_matrix['C Elemental'],
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
        
        organic_mass = (1.8 * concentration_matrix['C Orgánico'] +
                        0.6 * concentration_matrix['C Orgánico'].where(events['Event']=='SN', other=0) +
                        0.6 * concentration_matrix['C Orgánico'].where(events['Event']=='SP', other=0) +
                        0.6 * concentration_matrix['C Orgánico'].where(events['Event']=='S', other=0))
        uorganic_mass = (1.8 * uncertainty_matrix['C Orgánico'] +
                         0.6 * uncertainty_matrix['C Orgánico'].where(events['Event']=='SN', other=0) +
                         0.6 * uncertainty_matrix['C Orgánico'].where(events['Event']=='SP', other=0) +
                         0.6 * uncertainty_matrix['C Orgánico'].where(events['Event']=='S', other=0))
        
        elemental_C = concentration_matrix['C Elemental']
        uelemental_C = uncertainty_matrix['C Elemental']
        
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
                                    1.8 * uncertainty_matrix['C Orgánico'],
                                    uncertainty_matrix['C Elemental'],
                                    2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                    1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                    2.42 * uncertainty_matrix['Fe'] ], axis=0)
        
            
    if equation == 'Hand_2011':
        inorganic_ions = 1.375 * concentration_matrix['SO4'] + 1.29 * concentration_matrix['NO3']
        uinorganic_ions = np.linalg.norm( [ 1.375 * uncertainty_matrix['SO4'], 1.29 * uncertainty_matrix['NO3'] ] ,axis=0)
        
        organic_mass = (1.8 * concentration_matrix['C Orgánico'] +
                        0.6 * concentration_matrix['C Orgánico'].where(events['Event']=='SN', other=0) +
                        0.6 * concentration_matrix['C Orgánico'].where(events['Event']=='SP', other=0) +
                        0.6 * concentration_matrix['C Orgánico'].where(events['Event']=='S', other=0))
        uorganic_mass = (1.8 * uncertainty_matrix['C Orgánico'] +
                         0.6 * uncertainty_matrix['C Orgánico'].where(events['Event']=='SN', other=0) +
                         0.6 * uncertainty_matrix['C Orgánico'].where(events['Event']=='SP', other=0) +
                         0.6 * uncertainty_matrix['C Orgánico'].where(events['Event']=='S', other=0))
        elemental_C = concentration_matrix['C Elemental'] # Un solo elemento
        uelemental_C = uncertainty_matrix['C Elemental']
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
                                       1.8 * uncertainty_matrix['C Orgánico'],
                                       uncertainty_matrix['C Elemental'],
                                       2.2 * uncertainty_matrix['Al'], 2.49 * uncertainty_matrix['Si'],
                                       1.63 * uncertainty_matrix['Ca'], 1.94 * uncertainty_matrix['Ti'],
                                       2.42 * uncertainty_matrix['Fe'],
                                       1.8 * uncertainty_matrix['Cl'] ], axis=0)
        
        
    closure = sum(categories.values())
    categories['unexplained'] = concentration_matrix['PM2.5'] - closure
    ucategories['uunexplained'] = np.linalg.norm( [uclosure, uncertainty_matrix['PM2.5'] ], axis=0)
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