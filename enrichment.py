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

paper_year='1995'

if paper_year=='1964':
    enrichment = pd.read_csv(f'Taylor_1964.csv', skiprows=1, nrows=52, usecols=['Elementos', 'Crustal average']).dropna()
    enrichment.rename(columns={'Elementos':'Element', 'Crustal average':'Crustal average (ppmw)'}, inplace=True)
    enrichment['Crustal average (ppmw)'] = enrichment['Crustal average (ppmw)'].replace(',', '.', regex=True).astype('float')
    enrichment.dropna()
elif paper_year=='1995':
    enrichment = pd.read_csv('Taylor_1995.csv').dropna()


Ti_val_ref = enrichment['Crustal average (ppmw)'].loc[enrichment['Element']=='Ti']
enrichment['Ti factor'] = enrichment['Crustal average (ppmw)'] / Ti_val_ref.values

matrix_enriched = matrix.copy()
matrix_enriched.rename(columns={'Na total':'Na'}, inplace=True)


matrix_enriched = matrix_enriched[matrix_enriched.columns.intersection(list(enrichment['Element']))]
matrix_enriched = matrix_enriched.div(matrix_enriched['Ti'], axis=0)

print((enrichment['Ti factor'].loc[enrichment['Element']=='Na'].values))

for i in matrix_enriched.keys():
    matrix_enriched[i] = matrix_enriched[i]/(enrichment['Ti factor'].loc[enrichment['Element']==i].values[0])


matrix_enriched.to_csv('enrichment_factors.csv')
matrix_enriched.describe().to_csv('enrichment_factos_1964_statistics.csv')

for key in matrix_enriched.keys():
    fig, ax = plt.subplots()
    matrix_enriched[key].plot(ax=ax, style='.-', label=key)
    ax.axhline(1, linestyle='dashed', color='k')
    ax.set_yscale('log')
    ax.set_ylabel('EF')
    ax.legend()
    ax.set_title(f'Enrichment factor with respect to Ti for {key}')
    fig.savefig(f'images/enrichment_factors/taylor_{paper_year}/{key}_enrichment_factor.png')
    plt.close()

enrich_fact_desc = matrix_enriched.describe()
print(enrich_fact_desc)

fig, ax = plt.subplots()
ax.bar(enrich_fact_desc.keys(), height=enrich_fact_desc.loc['mean'], yerr=enrich_fact_desc.loc['std'])
ax.set_yscale('log')
ax.axhline(1, linestyle='dashed', color='k')
fig.savefig(f'images/enrichment_factors/taylor_{paper_year}/mean_enrichment.png')
plt.close()











