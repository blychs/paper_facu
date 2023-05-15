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
#     display_name: analysis
#     language: python
#     name: python3
# ---

# %%
##%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd

matrix = pd.read_excel('PMF_BA_full.xlsx', sheet_name="CONC", index_col=0)

matrix[matrix < 0] = np.nan
matrix = matrix.reset_index()
matrix["date"] = pd.to_datetime(matrix["date"])
matrix = matrix.rename(columns={'C Orgánico': 'OC',
                                 'C Elemental': 'EC',
                                 'C Total': 'TC',
                                 'Na sol': 'Na+',
                                 'Na no sol': 'InsolNa',
                                 'PM2,5': 'PM2.5'})
#print(matrix.keys())

keys = ['PM2.5', 'Cl', 'NO3', 'SO4', 'Na+', 'NH4', 'InsolNa',
       'Na total', 'Mg', 'Al', 'K', 'Ca', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
       'Ni', 'Cu', 'Zn', 'As', 'Se', 'Mo', 'Ag', 'Cd', 'Sb', 'Ba', 'Pb',
       'OCPk1 C', 'OCPk2 C', 'OCPk3 C', 'OCPk4 C', 'OCPk5 C', 'Pyrol C',
       'ECPk1 C', 'ECPk2 C', 'ECPk3 C', 'ECPk4 C', 'ECPk5 C', 'ECPk6 C', 'OC',
       'EC', 'TC']

clusters = pd.read_csv('Fechas Clusters.txt', delim_whitespace=True)
clusters = clusters.rename(
    columns={"YR": "year", "MO": "month", "DA": "day", "HR": "h"})

clusters["year"] = clusters["year"] + 2000

clusters["date"] = pd.to_datetime(clusters[["year", "day", "month", "h"]])

clusters = clusters[["date", "CL#"]].sort_values(by=["date"])
clusters["date"] = clusters["date"] - np.timedelta64(3, "h")
#print(matrix["date"].dt.date)
#print(clusters["date"].dt.date)
print(len(clusters))

clusters[keys] = 0

#print((matrix["date"].dt.date))

#print(clusters)
newmatrix = pd.DataFrame(np.repeat(matrix.values, 4, axis=0))
newmatrix.columns = matrix.columns

newmatrix["date"] = clusters["date"]
newmatrix["cluster"] = clusters["CL#"]

print(newmatrix.keys())
newmatrix = newmatrix.set_index(newmatrix["date"])
newmatrix = newmatrix.drop("date", axis=1)

# %%

add_noise = lambda x: x +  np.random.rand() * 0.2 - 0.1
for i in newmatrix.keys():
    fig, ax = plt.subplots()

    ax.scatter(newmatrix['cluster'].apply(add_noise), newmatrix[i],
                marker='.', c=newmatrix.index, cmap="tab20b")
    ax.set_xlabel("cluster")
    ax.set_ylabel(i)
    fig.savefig(f'images/cluster_elements/{i}.png')
    plt.close()

#%%
%matplotlib widget
#Explore elements
element = "Cl"
fig, ax = plt.subplots()
ax.scatter(newmatrix['cluster'].apply(add_noise), newmatrix[element],
            marker='.', c=newmatrix.index, cmap='tab20b')
ax.set_xlabel("Cluster")
ax.set_ylabel(element)
plt.show()