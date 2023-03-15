
def load_data(data_matrix, unc_matrix, gases, meteo, events):
    
    """
    Loads data for analysis of BA_2019
    """

    import datetime as dt
    import pandas as pd
    
    matrix = pd.read_excel(data_matrix, decimal=',', sheet_name='CONC')
    matrix = matrix.rename(columns={'PM2,5': 'PM2.5'})
    matrix['date'] = pd.to_datetime(matrix['date'])#.dt.date
    matrix.set_index(matrix['date'], inplace=True)
    matrix.drop('date', inplace=True, axis=1)
    matrix[matrix < 0] = np.nan
    matrix = matrix.reindex(sorted(matrix.columns), axis=1)
    
    unc = pd.read_excel(unc_matrix, decimal=',', sheet_name='UNC')
    unc = unc.rename(columns={'PM2,5': 'PM2.5'})
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
    events = pd.read_excel('BA_events.xlsx', index_col='date')
    
    gases = pd.read_csv('gases_mean.csv')
    gases['date'] = pd.to_datetime(gases['date'])
    gases.set_index(gases['date'], inplace=True)
    gases = matrix.join(gases)
    return matrix, unc, meteo, gases, events