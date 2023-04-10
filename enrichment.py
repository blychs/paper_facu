import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
from load_data import load_data

matrix, unc, meteo, gases, events = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                               'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                               'BA_events.xlsx')

enrichment = pd.read_csv('Taylor.csv', skiprows=1, nrows=52, usecols=['Elementos', 'Crustal average']).dropna()
enrichment.rename(columns={'Elementos':'Elements', 'Crustal average':'Crustal average (ppmw)'}, inplace=True)
enrichment.dropna()
enrichment['Crustal average (ppmw)'] = enrichment['Crustal average (ppmw)'].replace(',', '.', regex=True).astype('float')

Ti_val = enrichment['Crustal average (ppmw)'].loc[enrichment['Elements']=='Ti']
enrichment['Ti factor'] = enrichment['Crustal average (ppmw)'] / Ti_val.values


matrix_enriched = matrix.copy()
matrix_enriched.rename(columns={'Na total':'Na'}, inplace=True)
matrix_enriched = matrix_enriched[matrix.columns.intersection(list(enrichment['Elements']))]








