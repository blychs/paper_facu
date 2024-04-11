# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Python3 (analysis)
#     language: python
#     name: analysis
# ---

# %%
import pandas as pd
import numpy as np
import datetime as dt
from parserSFC import readSFC

def calc_OM(data_oc: pd.Series, event: pd.DataFrame, event_columname: str,
            event_labels: list[str], event_omoc: float | int = 2.6,
            noevent_omoc: float | int = 2.0) -> pd.Series:
    """
    Calculates the OM based on the OM/OC factor (beta_oc) for events and no events
    """
    data_noevent: pd.Series = noevent_omoc * data_oc.where(
        ~event[event_columname].isin(event_labels), other=0
    )
    data_event: pd.Series = event_omoc * data_oc.where(
        event[event_columname].isin(event_labels), other=0
    )
    return data_event + data_noevent


data: pd.DataFrame = pd.read_excel('PMF_BA_full.xlsx')
data['date'] = pd.to_datetime(data['date'])+ dt.timedelta(hours=12)
data = data.set_index(data['date'])
data = data.drop('date', axis=1)

events: pd.DataFrame = pd.read_excel('BA_events_testM.xlsx')
events['date'] = pd.to_datetime(events['date'] + dt.timedelta(hours=12))
events = events.set_index(events['date'])
events = events.drop('date', axis=1)
data = pd.concat([data, events[['Event']]], axis=1)
data["OM"] = calc_OM(data["C Orgánico"], events, event_columname="Event_M",
                     event_labels=["S", "SP", "SN","SL"])


meteo: pd. DataFrame = readSFC('OBSERVATORIO.SFC')

orig_hours: pd.Series = data.index


hour: int
for hour in orig_hours:
    for extended in range(1, 24):
        data.loc[hour + dt.timedelta(hours=extended)] = data.loc[hour]
data = data.sort_index()
meteo = meteo.loc[list(data.index)]

data = pd.concat([data, meteo[['ws', 'wd']]], axis=1)
print(data)
data[data.loc[:, data.columns!='Event'] < 0] = np.nan

data.to_csv('data_every_hour_obs_eventM.csv')

# %%
