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
import matplotlib as mpl
import matplotlib.cm as cm
import pandas as pd
import plotnine as pn
from statsmodels.multivariate.pca import PCA

matrix = pd.read_excel("PMF_BA_full.xlsx", sheet_name="CONC", index_col=0)

matrix[matrix < 0] = np.nan
matrix = matrix.reset_index()
matrix["date"] = pd.to_datetime(matrix["date"])
matrix = matrix.rename(
    columns={
        "C Orgánico": "OC",
        "C Elemental": "EC",
        "C Total": "TC",
        "Na sol": "Na$^+$",
        "Na no sol": "Na\nno sol",
        "PM2,5": "PM2.5",
    }
)
# print(matrix.keys())

keys = [
    "PM2.5",
    "Cl",
    "NO3",
    "SO4",
    "Na$^+$",
    "NH4",
    "Na\nno sol",
    "Na total",
    "Mg",
    "Al",
    "K",
    "Ca",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "As",
    "Se",
    "Mo",
    "Ag",
    "Cd",
    "Sb",
    "Ba",
    "Pb",
    "OCPk1 C",
    "OCPk2 C",
    "OCPk3 C",
    "OCPk4 C",
    "OCPk5 C",
    "Pyrol C",
    "ECPk1 C",
    "ECPk2 C",
    "ECPk3 C",
    "ECPk4 C",
    "ECPk5 C",
    "ECPk6 C",
    "OC",
    "EC",
    "TC",
]

clusters = pd.read_csv("Fechas Clusters.txt", delim_whitespace=True)
clusters = clusters.rename(
    columns={"YR": "year", "MO": "month", "DA": "day", "HR": "h"}
)

clusters["year"] = clusters["year"] + 2000

clusters["date"] = pd.to_datetime(clusters[["year", "day", "month", "h"]])

clusters = clusters[["date", "CL#"]].sort_values(by=["date"])
clusters["date"] = clusters["date"] - np.timedelta64(3, "h")
clusters = clusters.reset_index()
# with pd.option_context("display.max_rows", None):
#  print(clusters.date)
# print(matrix["date"].dt.date)
# print(clusters["date"].dt.date)

clusters[keys] = 0

# print((matrix["date"].dt.date))

# print(clusters)
newmatrix = pd.DataFrame(np.repeat(matrix.values, 4, axis=0))
newmatrix.columns = matrix.columns

newmatrix["date"] = clusters["date"]
newmatrix["cluster"] = clusters["CL#"]

newmatrix = newmatrix.set_index(newmatrix["date"])
newmatrix = newmatrix.drop("date", axis=1)

# %%


def add_noise(x):
    return x + np.random.rand() * 0.2 - 0.1


for i in newmatrix.keys():
    fig, ax = plt.subplots()

    ax.scatter(
        newmatrix["cluster"].apply(add_noise),
        newmatrix[i],
        marker=".",
        c=newmatrix.index,
        cmap="tab20b",
    )
    ax.set_xlabel("cluster")
    ax.set_ylabel(i)
    fig.savefig(f"images/cluster_elements/{i}.png")
    plt.close()

# %%
# %matplotlib widget
# Explore elements
element = "Cd"
fig, ax = plt.subplots()
ax.scatter(
    newmatrix["cluster"].apply(add_noise),
    newmatrix[element],
    marker=".",
    c=(newmatrix.index - np.timedelta64(12, "h")),
    cmap="tab20b",
)
ax.set_xlabel("Cluster")
ax.set_ylabel(element)
plt.show()
# %%
# %matplotlib widget
plt.plot(matrix["date"], matrix["InsolNa"])
# %%
l_todrop = [
    "cluster",
    "ECPk1 C",
    "ECPk2 C",
    "ECPk3 C",
    "ECPk4 C",
    "ECPk5 C",
    "ECPk6 C",
    "TC",
    "OCPk1 C",
    "OCPk2 C",
    "OCPk3 C",
    "OCPk4 C",
    "OCPk5 C",
    "Pyrol C",
    "PM2.5",
    "Na total",
]
newmatrix_1 = newmatrix[newmatrix["cluster"] == 1].drop(l_todrop, axis=1)
newmatrix_2 = newmatrix[newmatrix["cluster"] == 2].drop(l_todrop, axis=1)
newmatrix_3 = newmatrix[newmatrix["cluster"] == 3].drop(l_todrop, axis=1)
newmatrix_4 = newmatrix[newmatrix["cluster"] == 4].drop(l_todrop, axis=1)
newmatrix_5 = newmatrix[newmatrix["cluster"] == 5].drop(l_todrop, axis=1)


pca_c1 = PCA(newmatrix_1.dropna().T)
pca_c2 = PCA(newmatrix_2.dropna().T)
pca_c3 = PCA(newmatrix_3.dropna().T)
pca_c4 = PCA(newmatrix_4.dropna().T)
pca_c5 = PCA(newmatrix_5.dropna().T)


pca_c1.plot_scree(log_scale=False)

fig, ax = plt.subplots()  # figsize=(5, 5))
pca_c1.loadings.plot.scatter(x="comp_00", y="comp_01", ax=ax, label="C#1", marker="1")
pca_c2.loadings.plot.scatter(
    x="comp_00", y="comp_01", color="k", ax=ax, label="C#2", marker="x"
)
pca_c3.loadings.plot.scatter(
    x="comp_00", y="comp_01", color="r", ax=ax, label="C#3", marker="+"
)
pca_c4.loadings.plot.scatter(
    x="comp_00", y="comp_01", color="g", ax=ax, label="C#4", marker="*"
)
pca_c5.loadings.plot.scatter(
    x="comp_00", y="comp_01", color="m", ax=ax, label="C#5", marker="2"
)
ax.set_xlabel("PCA 1")
ax.set_ylabel("PCA 2")
plt.plot()

# %%
compare_PCA = pd.concat(
    [
        pca_c1.factors["comp_00"],
        pca_c2.factors["comp_00"],
        pca_c3.factors["comp_00"],
        pca_c4.factors["comp_00"],
        pca_c5.factors["comp_00"],
    ],
    axis=1,
)

# display(compare_PCA)

for i in range(0, 3):
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(pca_c1.factors.iloc[:, i], lw=1, alpha=0.6, label="PC C1")
    ax.plot(pca_c2.factors.iloc[:, i], lw=1, alpha=0.6, label="PC C2")
    ax.plot(pca_c3.factors.iloc[:, i], lw=1, alpha=0.6, label="PC C3")
    ax.plot(pca_c4.factors.iloc[:, i], lw=1, alpha=0.6, label="PC C4")
    ax.plot(pca_c5.factors.iloc[:, i], lw=1, alpha=0.6, label="PC C5")
    ax.axhline(0, color="k", linestyle="dotted")
    ax.set_xticklabels(newmatrix_1.columns.values, size=10)
    # ax.set_xlim(0, 51)
    ax.set_xlabel("Species", size=17)
    ax.set_ylabel(f"PCA{i+1} composition")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"images/cluster_elements/pca{i+1}_comp.png")
    plt.close()
# %%
pca_all = PCA(newmatrix.drop(l_todrop, axis=1).dropna().T)
my_cmap = mpl.colors.LinearSegmentedColormap.from_list(
    "my_own_cmap", [(0, 1, 0.5), (1, 0.5, 0), (0, 0, 1), (0, 0, 0), (1, 0, 1)], N=5
)
# display(pca_all.factors)
pca_all_withcluster = pd.concat([pca_all.loadings, newmatrix["cluster"]], axis=1)
pca_all_withcluster.plot.scatter(
    x="comp_00", y="comp_01", c=pca_all_withcluster["cluster"], cmap=my_cmap
)  # ,
plt.xlabel("PC 1")
plt.ylabel("PC 2")
plt.show()

# %%
# display(newmatrix)
fig, ax = plt.subplots(figsize=(10, 5))
print(newmatrix.index)

ax.plot(newmatrix.index, newmatrix["Cu"], ":k")
ax.plot(
    newmatrix.index, newmatrix["Cu"].where(newmatrix["cluster"] == 1), ".-", label="C1"
)
ax.plot(
    newmatrix.index, newmatrix["Cu"].where(newmatrix["cluster"] == 2), ".-", label="C2"
)
ax.plot(
    newmatrix.index, newmatrix["Cu"].where(newmatrix["cluster"] == 3), ".-", label="C3"
)
ax.plot(
    newmatrix.index, newmatrix["Cu"].where(newmatrix["cluster"] == 4), ".-", label="C4"
)
ax.plot(
    newmatrix.index, newmatrix["Cu"].where(newmatrix["cluster"] == 5), ".-", label="C5"
)
ax.legend()
ax.set_xlabel("Date")
ax.set_ylabel("Cu (µg/m$^3$)")
fig.tight_layout()
plt.show()


# %%
