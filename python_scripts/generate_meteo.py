# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import pandas as pd
import datetime as dt

# +
#Carga meteorologia
obsbaires = pd.read_csv('../../doctorado/paper_laura/ARCAL_ISH/obsbaires.csv', delimiter=';') # obsbaires2 es abrir y volver a cerrar en excel para que ponga bien los delimiters, si no mezcla ; y ,
obsbaires = obsbaires[(obsbaires['date[yyyymmddHHMM]'] >= 201904031200)]
fechas = df['Date']
# Reemplazo los valores que son nan
obsbaires['wdir'] = obsbaires['wdir'].where(obsbaires['wdir'] < 999)
obsbaires['wspd[m/s]'] = obsbaires['wspd[m/s]'].where(obsbaires['wspd[m/s]'] < 999)
obsbaires['clht[km]'] = obsbaires['clht[km]'].where(obsbaires['clht[km]'] < 99)
obsbaires['dptp[C]'] = obsbaires['dptp[C]'].where(obsbaires['dptp[C]'] < 999)
obsbaires['slvp[hPa]'] = obsbaires['slvp[hPa]'].where(obsbaires['slvp[hPa]'] < 9999)
obsbaires['press[hPa]'] = obsbaires['press[hPa]'].where(obsbaires['press[hPa]'] < 9999)
obsbaires['prcp[mm]'] = obsbaires['prcp[mm]'].where(obsbaires['prcp[mm]'] < 999)
obsbaires['sky[octas]'] = obsbaires['sky[octas]'].where(obsbaires['sky[octas]'] > 99)
obsbaires['skyOpaque[octas]'] = obsbaires['skyOpaque[octas]'].where(obsbaires['skyOpaque[octas]'] < 99)
obsbaires['tmpd[C]'] = obsbaires['tmpd[C]'].where(obsbaires['tmpd[C]'] < 99)

# Calculo RH como rh=100*(EXP((17.625*d)/(243.04+d))/EXP((17.625*t)/(243.04+t)));
# donde d=dewpoint (dptp) y t=dry temp (tmpd), aprox de August-Roche-Magnus
d = obsbaires['dptp[C]']
t = obsbaires['tmpd[C]']
obsbaires['RH[%]'] = 100 * (np.exp((17.625 * d)/(243.04 + d)) / np.exp((17.625 * t)/(243.04 + t)))
# 

#Promedio fechas que están dobles
obsbaires = obsbaires.groupby('date[yyyymmddHHMM]').mean().reset_index()



# Parsing de fecha
datetime = pd.to_datetime(obsbaires['date[yyyymmddHHMM]'], format='%Y%m%d%H%M')
obsbaires['date'] = datetime
# Reordeno columnas-
cols = obsbaires.columns.tolist()
cols = cols[-1:] + cols[:-1] # Mando la última columna (date) al principio
obsbaires = obsbaires[cols]
#display(obsbaires)
obsbaires = obsbaires.set_index('date')

#display(obsbaires)

obsbaires = obsbaires.groupby(level=0).sum()
#display(obsbaires)

obsbaires24hs_mean = obsbaires.resample('24H', offset='12h').mean() # Resampleo para tener promedios diarios

#display(obsbaires24hs_mean)
obsbaires24hs_sum = obsbaires.resample('24H', offset='12h').sum() # Resampleo para tener promedios diarios

obsbaires_comb = pd.concat([obsbaires24hs_mean[['wdir', 'wspd[m/s]', 'clht[km]', 'hzvs[km]', 'tmpd[C]',
                                          'slvp[hPa]', 'press[hPa]', 'sky[octas]', 'skyOpaque[octas]', 'RH[%]']],
                         obsbaires24hs_sum[['prcp[mm]', 'prcpPeriod[hours]']]], axis=1)

#display(obsbaires_comb)
obsbaires_comb.index = obsbaires_comb.index.normalize()
obsbaires_comb = (obsbaires_comb.reindex(index = fechas.to_list()))
obsbaires_comb.to_excel('obsbaires_meteo2.xlsx')
# 
# #with pd.
# #display(fechas['obsbaires'])
