# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
# #%matplotlib widget

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
from parserSFC import readSFC, ventilation_coef

df = pd.read_csv('../paper_laura/ARCAL_ISH/obsbaires.csv', delimiter=';')
df = df.rename({'date[yyyymmddHHMM]':'date', 'wspd[m/s]':'ws', 'wdir':'wd'}, axis=1)
df['date'] = pd.to_datetime(df['date'].astype(str))
df = df[df['date'] > dt.datetime(2019,4,3)]
df = df[df['date'] < dt.datetime(2020,4,1)]
df_filters = pd.read_excel('PMF_BA_full.xlsx', sheet_name='CONC')
df_filters['PM2,5'] = df_filters['PM2,5'].where(df_filters['PM2,5'] > 0)
print((df['ws'] == 0).sum()/len(df['ws']) * 100)
# Index and select correct timestamps
df = (df.loc[(df['date'].dt.date.isin(df_filters['date'].dt.date) & (df['date'].dt.hour >= 12)) |
              (df['date'].dt.date.isin(df_filters['date'].dt.date + dt.timedelta(days=1)) & (df['date'].dt.hour < 12))])

events = pd.read_excel('BA_events.xlsx', index_col='date')
display(events)

#with pd.option_context('display.max_rows', None,):
#    display(df)


    
def num_calms_samplingdate(df, offset=-12):
    """ Calculate the number of calm
    hours per sampling date, with and offset of
    a number of hours"""
    df2 = df.filter(['date', 'ws'], axis=1) # copy only needed columns
    df2['date'] = df2['date'] + dt.timedelta(hours=offset)
    df2['calms'] = (df['ws'] == 0.0)
    num_of_calms = df2['calms'].groupby(df2['date'].dt.date).sum()
    return(num_of_calms)

def average_windspeed(df, offset=-12):
    """ Calculate the average windspeed,
    offset of a number of hours"""
    df2 = df.filter(['date', 'ws'], axis=1) # copy only needed columns
    df2['date'] = df2['date'] + dt.timedelta(hours=offset)
    ws = df2['ws'].groupby(df2['date'].dt.date).mean()
    return(ws)


num_calms = num_calms_samplingdate(df)

    
#fig, ax = plt.subplots()
#ax.plot(df['date'], df['ws'], '.')
#ax.set_xlabel('Date')
#ax.set_ylabel('Wind speed (m/s)')
#plt.show()
daily_ws = average_windspeed(df)

#print(daily_ws)
#print(num_calms)

data_out = pd.DataFrame(daily_ws)

data_out['num_calms'] = num_calms
display(data_out)
data_out.to_excel('wind_data.xlsx')
#data_out = data_out.set_index(daily_ws.index)

#with plt.style.context('ggplot'):

# #%matplotlib widget
fig, ax = plt.subplots()
#num_calms.plot(ax=ax, style='.-')
#ax.axhline(6.58, color='k', linestyle='dashed', label='average')
ax.plot(df_filters['date'], df_filters['PM2,5'], '.-')
ax.bar(num_calms.index, num_calms * 5, color='k', label='Number of calm hours')
#ax.plot(events.index, )
ax.set_ylabel('Number of hours registered as calm')
plt.show()

fig, ax = plt.subplots()
daily_ws.plot(ax=ax, style='.-', label='Wind speed', zorder=1)
# ax.plot(df['date'], df['ws'], zorder=0)
ax.bar(num_calms.index, num_calms, color='k', label='Number of calm hours')
#ax.axhline(6.58, color='k', linestyle='dashed', label='average')
ax.set_ylabel('Average ws (m/s)')
#plt.xticks(rotation=90)
plt.legend()
plt.show()
    


# + tags=[]
df_uaest = readSFC('OBSERVATORIO.SFC')
df_uaest = (df_uaest.loc[(df_uaest['date'].dt.date.isin(df_filters['date'].dt.date) & (df_uaest['date'].dt.hour.astype(int) >= 12 )) |
              (df_uaest['date'].dt.date.isin(df_filters['date'].dt.date + dt.timedelta(days=1)) & (df_uaest['date'].dt.hour.astype(int) < 12))])

#display(df_uaest)
#with pd.option_context('display.max_rows', None):#, 'display.max_columns', None,):
#    display(df_uaest[['ws']])
#display(df[['ws']])

# +
df_uaest['VentCoef'] = ventilation_coef(df_uaest)
with pd.option_context('display.max_columns', None, 'display.max_rows', None):
    display(df_uaest[['ws', 'PBLHc', 'PBLHm', 'L', 'VentCoef']])

df_uaest.to_csv('datos_meteo_obs_filtro.csv')
