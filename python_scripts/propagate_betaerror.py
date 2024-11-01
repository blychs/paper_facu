#%%
import numpy as np
import pandas as pd
from funciones_pmfBA import estimation_om_oc
from load_data import load_data

def propag_error_in_average(error_list, events=None, event_label=None):
    """
    Progpagate the errors in beta
    """

    return np.linalg.norm(np.array(error_list), ord=1)/len(error_list)



matrix, unc, meteo, gases, events = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                              'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                              'BA_events_testM.xlsx')

methods = ['Macias_1981', 'Solomon_1989', 'Chow_1994',
           'Malm_1994', 'Chow_1996', 'Andrews_2000',
           'Malm_2000', 'Maenhaut_2002', 'DeBell_2006',
           'Hand_2011','Simon_2011']

event_columnname="Event_M"
event_labels= ["S", "SP", "SN","SL"]

l_errors_all = []
l_errors_events = []
l_errors_no_events = []

for method in methods:
    om_oc = estimation_om_oc(matrix, method=method)
    om_oc_events = estimation_om_oc(
        matrix.where(events[event_columnname].isin(event_labels)),
        method=method
        )
    om_oc_no_events = estimation_om_oc(
        matrix.where(~events[event_columnname].isin(event_labels)),
        method=method
        )
    l_errors_all.append(om_oc.bse[1])
    l_errors_events.append(om_oc_events.bse[1])
    l_errors_no_events.append(om_oc_no_events.bse[1])

print(l_errors_all)
print(l_errors_events)
print(l_errors_no_events)
print(f"error OM/OC all = {propag_error_in_average(l_errors_all):.1g}")
print(f"error OM/OC event = {propag_error_in_average(l_errors_events):.1g}")
print(f"error OM/OC no event = {propag_error_in_average(l_errors_no_events):.1g}")