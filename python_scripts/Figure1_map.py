#%%
directory='/home/usuario/mdiaz/Documents/paper_facu/python_scripts'
os.chdir(directory)

import cartopy
import cartopy.crs as ccrs
from cartopy.io.img_tiles import OSM
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patheffects
from matplotlib.ticker import StrMethodFormatter
import numpy as np
import datetime
import geopandas as gpd
from module_maps import scale_bar
import os 


directory='/home/usuario/mdiaz/Documents/paper_facu'
os.chdir(directory)
proyecto="CNEA"
epsg="32721"
utmzone=21


predioPath=directory+"/GIS/predioCNEA.csv"
predio = gpd.read_file(predioPath)
predio.crs = 'epsg:'+epsg

dominioPath=directory+"/GIS/dominioCNEA.csv"
dominio = gpd.read_file(dominioPath)
dominio.crs = 'epsg:'+epsg

cneaPath=directory+"/GIS/puntoCNEA.csv"
cnea = gpd.read_file(cneaPath)
cnea.crs = 'epsg:'+epsg

edificio48Path=directory+"/GIS/edificio48v3.csv"
edificio48 = gpd.read_file(edificio48Path)
edificio48.crs = 'epsg:'+epsg

#request=cimgt.GoogleTiles()#style='satellite')
projection = ccrs.UTM(zone=utmzone,southern_hemisphere=True)
request = OSM()

BORDERS2_10m = cartopy.feature.NaturalEarthFeature('cultural', 'admin_1_states_provinces','10m', edgecolor='black', facecolor='none')

def blank_axes(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')
    ax.tick_params(labelbottom='off', labeltop='off', labelleft='off', labelright='off',bottom='off', top='off', left='off', right='off' )
    ax.axis('off')
#end blank_axes

fig = plt.figure(figsize=(6.4,4.8))#figsize=(10,12))
# ------------------------------- Surrounding frame ------------------------------
# set up frame full height, full width of figure, this must be called first
left = 0.0
bottom = -0.05
width = 1.0
height =1.0
rect = [left,bottom,width,height]
ax4 = plt.axes(rect)

# turn on the spines we want, ie just the surrounding frame
blank_axes(ax4)
ax4.spines['right'].set_visible(True)
ax4.spines['top'].set_visible(True)
ax4.spines['bottom'].set_visible(False)
ax4.spines['left'].set_visible(False)

#ax4.text(0.01,0.01,' Sistema espacial de referencia (SRC) del proyecto EPSG:'+epsg+'. Mapa generado el '+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+' para '+proyecto, fontsize=8)
#ax4.text(0.01,0.01,' EPSG:'+epsg+'. ')#'. Mapa generado el '+datetime.datetime.now().strftime("%Y-%m-%d")+' para '+proyecto, fontsize=8)

# -------------------------Mapa Predio (Main)----------------------
# set up main map almost full height (allow room for title), right 80% of figure
left = 0.323
bottom = -0.002
width = 0.76
height = 0.905
rect = [left,bottom,width,height]
#ax1 = plt.axes(rect, projection=ccrs.PlateCarree(), )
ax1 = plt.axes(rect, projection=projection)

lon0,lat0,lon1,lat1=dominio.geometry.total_bounds
borde=3700.0
borde2=600.0
ax1.set_extent((lon0+borde+borde2,lon1-borde-borde2,lat0+borde,lat1-borde), crs=projection)
ax1.add_image(request,14)#, interpolation='spine36')
#ax1.add_image(request,8)#, interpolation='spine36')

predio.plot(ax=ax1,color="none",edgecolor="red",alpha=0.9,linewidth=1.1)
dominio.plot(ax=ax1, color="none",edgecolor="green", alpha=0.9,linewidth=1.1)
edificio48.plot(ax=ax1, color="none",edgecolor="blue", alpha=0.9,linewidth=1.1)

#cnea.plot(ax=ax1, color="blue",marker="o", alpha=0.9,linewidth=1.1)
#cnealatlon=cnea.to_crs(4326)
#x = cnealatlon.centroid.x
#y = cnealatlon.centroid.y
plt.scatter(361250,6173280,color="green",alpha=0.9)
#plot(x,y, color='blue', marker='o', linestyle='dashed', linewidth=2, markersize=12)
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top')
ax1.yaxis.tick_right()
ax1.yaxis.set_label_position("right")
#ax1.set_xticks(np.linspace(lon0+borde+borde2,lon1-borde-borde2,3).round(0), crs=projection)
#ax1.set_yticks(np.linspace(lat0+borde,lat1-borde,4).round(0) , crs=projection)
#ax1.tick_params(labelsize=8)
#fig.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}')) # No decimal places on axis


#### bar offset is how far from bottom left corner scale bar is (x,y) and how far up is scale bar text
bar_offset = [0.1, 0.01, 0.1]
#bar_offset = [0.0, 0.0, 0.0]
bar_lon0 = lon0+borde+borde2 + (lon1-lon0)*bar_offset[0]
bar_lat0 = lat0+borde + (lat1-lat0)*bar_offset[1]
#
text_lon0 = bar_lon0
text_lat0 = lat0+borde + (lat1-lat0)*bar_offset[2]
bar_tickmark = 1000 # metres
bar_ticks = 5
bar_alpha = 0.7#0.3
#
bar_color = ["black","white"]#['black', 'red']
#
buffer = [patheffects.withStroke(linewidth=3, foreground="w")]

### Plot the scalebar label
units = 'm'
#t0 = ax1.text(text_lon0, text_lat0, str(bar_tickmark/1000) + ' ' + units, horizontalalignment='left', verticalalignment='bottom', path_effects=buffer, zorder=2)
scale_bar(ax1, 500)

# draw a scale bar that is a set of colored line segments (bar_ticks of these), bar_tickmarks long
for i in range(bar_ticks):
    end_lat, end_lon = bar_lat0, bar_lon0+bar_tickmark/bar_ticks
    #TODO make transform match ax projection
    ax1.plot([bar_lon0, end_lon], [bar_lat0, end_lat], color=bar_color[i%2], linewidth=7, solid_capstyle='butt', alpha = bar_alpha)
    # start of next bar is end of last bar
    bar_lon0 = end_lon
    bar_lat0 = end_lat
#end for
#
# ---------------------------------Mapa Dominio------------------------
left = 0
bottom = 0
width =  0.32
height = 0.47
rect = [left,bottom,width,height]




#ax2 = plt.axes(rect, projection=ccrs.PlateCarree(), )
ax2 = plt.axes(rect, projection=projection )

lon0,lat0,lon1,lat1=dominio.geometry.total_bounds
borde=29000
ax2.set_extent((lon0-borde,lon1+borde,lat0-borde,lat1+borde), crs=projection)

ax2.add_image(request,11)#, interpolation='spine36')
#ax2.add_image(request,5)#, interpolation='spine36')

dominio.plot(ax=ax2,color="none",edgecolor="green",alpha=0.9)
predio.plot(ax=ax2,color="none",edgecolor="red",alpha=0.9)
cnea.plot(ax=ax2,color="blue")
## --------------------------------Mapa Argentina------------------------
left = 0
bottom = 0.42
width =  0.32
height = 0.495
rect = [left,bottom,width,height]

ax3 = plt.axes(rect, projection=ccrs.PlateCarree() )
#ax3.set_extent((-78,-52,-57,-10))
ax3.set_extent((-83,-42,-55,-15))
ax3.stock_img()
ax3.add_feature(BORDERS2_10m, edgecolor='grey', linewidth=0.8)
#blank_axes(ax3)

domLatLon=dominio.to_crs(4326)
x = domLatLon.centroid.x
y = domLatLon.centroid.y
plt.scatter(x,y,color="green",alpha=0.9)#,s=1.7,marker="o")
## ---------------------------------North Arrow  ----------------------------
left = 0.95
bottom = 0.78
width = 0.2
height = 0.2
rect = [left,bottom,width,height]
ax5 = plt.axes(rect)

# need a font that support enough Unicode to draw up arrow. need space after Unicode to allow wide char to be drawm?
ax5.text(0.5, 0.0,u'\u25B2 \nN ', ha='center', fontsize=18, rotation = 0,alpha=0.7)
blank_axes(ax5)
#plt.tight_layout()
#plt.show()

## --------------------------Save figure -----------------------------
plt.savefig(directory+'/images/plotUbicacionv4.png',dpi=300, bbox_inches = 'tight',pad_inches = 0.1)
plt.close()
# %%
