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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

data = pd.read_csv('CNEA-CAC-AQdata20192020.csv')
data["date"]=pd.to_datetime(data["date"])
data.set_index(data["date"], inplace=True)

df_filters = pd.read_excel('PMF_BA_full.xlsx', sheet_name='CONC')
df_filters['PM2,5'] = df_filters['PM2,5'].where(df_filters['PM2,5'] > 0)

data = (data.loc[(data['date'].dt.date.isin(df_filters['date'].dt.date) & (data['date'].dt.hour.astype(int) >= 12 )) |
              (data['date'].dt.date.isin(df_filters['date'].dt.date + dt.timedelta(days=1)) & (data['date'].dt.hour.astype(int) < 12))])
data_daily = data.resample('24H', offset='12H').mean(numeric_only=True).dropna()

data_daily.index = data_daily.index.to_period('D')
data_daily[['NO', 'NO2', 'NOx', 'CO', 'O3', 'SO2']].to_csv('gases_mean.csv')
#display(data)
