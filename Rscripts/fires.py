import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime, timedelta

# Configuración de rutas
pathgraphs = "Figures"
upper_windspeed = 8.5
pscflims = (0, 1)
colorset = "plasma"
pathtraj = "../tdumps"
pattern = "tdump*"
bbpdatafile = "../data/data_every_hour_obsv5.csv"

# Cargar dataset de calidad del aire
bbpdata = pd.read_csv(bbpdatafile)
bbpdata['date'] = pd.to_datetime(bbpdata['date'], utc=False).dt.tz_localize(None)
import os
import re
import pandas as pd
from datetime import datetime, timedelta

def read_tdumps(hours=72, pathtraj="", height=None, location=None, pattern="tdump*", 
                pathtrajout="", onlyseason=None, outfilename="trajcombined.txt", numheaderlines=13):
    """
    Lee archivos tdump de HYSPLIT, filtra por parámetros específicos y los combina en un DataFrame.
    """
    # Obtener lista de archivos según el patrón
    files = [os.path.join(pathtraj, f) for f in os.listdir(pathtraj) if re.search(pattern, f)]
    
    # Filtrar por ubicación si se especifica
    if location:
        files = [f for f in files if location in f]
    
    # Filtrar por temporada si se especifica
    if onlyseason:
        datepattern = re.compile(r"(\d{2})(\d{4})" + location)
        filtered_files = []
        for f in files:
            match = datepattern.search(f)
            if match:
                month = int(match.group(1))
                seasonnames = ['V', 'O', 'I', 'P']
                indx = {12: 'V', 1: 'V', 2: 'V', 3: 'O', 4: 'O', 5: 'O',
                        6: 'I', 7: 'I', 8: 'I', 9: 'P', 10: 'P', 11: 'P'}
                if indx[month] == onlyseason:
                    filtered_files.append(f)
        files = filtered_files

    combined_data = []
    for file in files:
        with open(file, 'r') as f:
            lines = f.readlines()
        
        # Encontrar el número de líneas del encabezado
        for i, line in enumerate(lines):
            if "     1 PRESSURE" in line:
                numheaderlines = i
                break
        
        # Leer datos eliminando el encabezado
        data = [line.split() for line in lines[numheaderlines+1:]]
        df = pd.DataFrame(data)
        combined_data.append(df)

    # Concatenar todos los datos en un solo DataFrame
    traj = pd.concat(combined_data, ignore_index=True)
    
    # Renombrar columnas
    traj.columns = ['receptor', 'drop1', 'year', 'month', 'day', 'hour', 
                    'drop2', 'drop3', 'hour_inc', 'lat', 'lon', 'height', 'pressure']
    traj = traj.drop(columns=['drop1', 'drop2', 'drop3'])

    # Ajustar el formato del año
    traj['year'] = traj['year'].astype(int)
    traj['year'] = traj['year'].apply(lambda y: y + 2000 if y < 50 else y + 1900)

    # Crear fechas y ajustar el tiempo de llegada
    traj['date2'] = pd.to_datetime(traj[['year', 'month', 'day', 'hour']])
    traj['date'] = traj['date2'] - pd.to_timedelta(traj['hour_inc'].astype(int), unit='h')

    # Crear una columna de fecha para fusionar
    traj['datemerge'] = traj['date'].dt.strftime('%d%m%Y')
    traj['datemerge'] = pd.to_datetime(traj['datemerge'], format='%d%m%Y')

    return traj

def add_chem(pollutant, traj):
    """
    Fusiona los datos de trayectoria con los datos de contaminación química.
    """
    traj['datemerge'] = traj['date'].dt.strftime('%d%m%Y')
    traj['datemerge'] = pd.to_datetime(traj['datemerge'], format='%d%m%Y')
    
    trajconchem = traj.merge(pollutant, on='datemerge', how='left')
    return trajconchem


traj500 = read_tdumps(pathtraj = pathtraj, pattern=pattern)
# bbpdata['date'] = pd.to_datetime(bbpdata['date'], utc=False).dt.tz_localize(None)
# traj500['puestadefiltro'] = traj500['date'] - pd.Timedelta(seconds=123600)
# traj500 = traj500[traj500['receptor'] == 3]

# # Fusionar trayectorias con datos químicos
# trajconchem = traj500.merge(bbpdata, on='date', how='left')
# rename_dict = {'C.Elemental': 'EC', 'C.Orgánico': 'OC', 'C.Total': 'TC', 'Na.sol': 'Na'}
# trajconchem = trajconchem.rename(columns=rename_dict)

# # Cargar datos de incendios
# fires = pd.read_csv("../data/fires/DL_FIRE_M-C61_563912/fire_archive_M-C61_563912.csv")
# fires['acq_time'] = fires['acq_time'].apply(lambda x: f"{int(x):04d}")
# fires['hour'] = fires['acq_time'].str[:2]
# fires['minute'] = fires['acq_time'].str[2:]
# fires['datetime'] = pd.to_datetime(fires['acq_date'] + ' ' + fires['hour'] + ':' + fires['minute'])

# # Filtrar incendios dentro del rango de fechas
# start_date = datetime(2019, 8, 26)
# end_date = datetime(2019, 8, 30)
# start_filter = start_date - timedelta(hours=72)
# end_filter = end_date
# fires_filtered = fires[(fires['datetime'] >= start_filter) & (fires['datetime'] <= end_filter)]

# # Graficar trayectorias y puntos de incendios
# fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': ccrs.SouthPolarStereo()})
# ax.add_feature(cfeature.LAND, facecolor='lightgray')
# ax.add_feature(cfeature.BORDERS, linestyle=':')
# ax.add_feature(cfeature.COASTLINE)
# ax.set_extent([-70, -50, -40, -20], crs=ccrs.PlateCarree())

# # Agrupar por fecha y trazar trayectorias con líneas
# for date, group in trajconchem.groupby('date'):
#     ax.plot(group['lon'], group['lat'], color='blue', transform=ccrs.PlateCarree(), alpha=0.5, linewidth=1)

# ax.scatter(fires_filtered['longitude'], fires_filtered['latitude'], s=10, color='red', transform=ccrs.PlateCarree(), label='Fires')
# ax.set_title("Trayectorias y Fuegos")
# ax.legend()
# plt.show()
