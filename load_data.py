
def load_data(data_matrix, unc_matrix, gases, meteo, events, clusters=None):
    
    """
    Loads data for analysis of BA_2019
    """

  
    import pandas as pd
    import numpy as np
    
    matrix = pd.read_excel(data_matrix, decimal=',', sheet_name='CONC')
    matrix = matrix.rename(columns={'PM2,5': 'PM2.5', 'C Orgánico': 'OC',
                                     'C Elemental': 'EC', 'C Total': 'TC'})
    matrix['date'] = pd.to_datetime(matrix['date'])#.dt.date
    matrix.set_index(matrix['date'], inplace=True)
    matrix.drop('date', inplace=True, axis=1)
    matrix[matrix < 0] = np.nan
    matrix = matrix.reindex(sorted(matrix.columns), axis=1)
    
    unc = pd.read_excel(unc_matrix, decimal=',', sheet_name='UNC')
    unc = unc.rename(columns={'PM2,5': 'PM2.5', 'C Orgánico': 'OC',
                                     'C Elemental': 'EC', 'C Total': 'TC'})
    unc['date'] = pd.to_datetime(unc['date'])#.dt.date
    unc.set_index(unc['date'], inplace=True)
    unc.drop('date', inplace=True, axis=1)
    unc[unc < 0] = np.nan
    unc = unc.reindex(sorted(unc.columns), axis=1)
    
    meteo = pd.read_csv(meteo)
    meteo = meteo.rename(columns={'Unnamed: 0': 'date'})
    meteo['date'] = pd.to_datetime(meteo['date'])
    meteo.set_index(meteo['date'], inplace=True)
    meteo.drop('date', inplace=True, axis=1)
    
    matrix = matrix.join(meteo)
    matrix['temp'] = (matrix['temp']-273.15).round(2)
    events = pd.read_excel(events, index_col='date')
    
    gases = pd.read_csv('gases_mean.csv')
    gases['date'] = pd.to_datetime(gases['date'])
    gases.set_index(gases['date'], inplace=True)
    gases = matrix.join(gases)
    if clusters == None:
        return matrix, unc, meteo, gases, events
    else:
        clusters = pd.read_csv(clusters)
        clusters["date"] = pd.to_datetime(clusters["date"])
        clusters.set_index(clusters["date"], inplace=True)
        clusters.drop("date", inplace=True, axis=1)
        return matrix, unc, meteo, gases, events, clusters