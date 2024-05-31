import pandas as pd
import numpy as np
import datetime as dt


data: pd.DataFrame = pd.read_excel("PMF_BA_fullv3.xlsx")
data["date"] = pd.to_datetime(data["date"]) + dt.timedelta(hours=12)
data = data.set_index(data["date"])
data = data.drop("date", axis=1)

events: pd.DataFrame = pd.read_excel("BA_events_testM.xlsx")
events["date"] = pd.to_datetime(events["date"] + dt.timedelta(hours=12))
events = events.set_index(events["date"])
events = events.drop("date", axis=1)
data = pd.concat([data, events[["Event_M"]]], axis=1)
data = data.drop(columns=["Na total", "Na no sol"])
data = data.rename(columns={"Na sol": "Na", "PM2,5": "PM2.5"})


meteo: pd.DataFrame = pd.read_csv("datos_meteo_blhera5_hourly.csv")
meteo["date"] = pd.to_datetime(meteo["date"])
print(meteo)
print(meteo.dtypes)

orig_hours: pd.Series = data.index


hour: int
for hour in orig_hours:
    for extended in range(1, 24):
        data.loc[hour + dt.timedelta(hours=extended)] = data.loc[hour]
data = data.sort_index()
data = data.reset_index()
print(data.dtypes)

data: pd.DataFrame = pd.merge(data, meteo, on="date")
print(data)
# data[data.loc[:, data.columns != "Event"] < 0] = np.nan

data.to_csv("data_every_hour_obs_eventM.csv", index=False, float_format="%.4g")

# %%
