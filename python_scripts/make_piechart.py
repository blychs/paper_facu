"""
HORRIBLE, HORRIBLE CODE
"""

import matplotlib as mpl
import cmocean
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import colorcet as cc

PATH = "../data"


def main() -> None:
    data: pd.DataFrame = pd.read_excel(
        f"{PATH}/PMF.xlsx",
        index_col=0,
        header=0,
        usecols=list(range(0, 8)),
        skiprows=[1],
    )

    data2: pd.DataFrame = pd.read_excel(
        f"{PATH}/72g_Constrained.xlsx",
        sheet_name="Profiles",
        skiprows=39,
        usecols=range(1, 9),
        nrows=33,
        header=None,
    )
    [data2.rename(columns={x + 1: f"Factor {x}"}, inplace=True)
     for x in range(1, 8)]
    data2 = data2.rename(columns={1: "Species"})
    data2 = data2[
        [
            "Species",
            "Factor 3",
            "Factor 6",
            "Factor 1",
            "Factor 7",
            "Factor 4",
            "Factor 5",
            "Factor 2",
        ]
    ]
    data2 = data2.set_index("Species", drop=True)
    data2 = data2.rename(
        index={"PM2,5": "PM$_{2.5}$", "Na sol": "Na", "Pyrol C": "Pyrol"}
    )
    for i in data2.index:
        if "Pk" in i:
            data2 = data2.rename(
                index={i: i.replace("Pk", "").replace(" C", "")})
    print(data2.dtypes)

    new_data = data.T.reset_index(level=0).rename(columns={"index": "Source"})
    new_data = new_data.sort_values(
        by="PM2,5", ascending=False, ignore_index=True)
    new_data["number"] = list(range(1, len(new_data["Source"]) + 1))
    new_data["number"] = new_data["number"].astype("str")
    new_data["Source"] = (new_data["Source"] + " ")
    new_data.loc[0, "Source"] = new_data["Source"][0].replace(" Bu", "\nBu")
    new_data.loc[1, "Source"] = new_data["Source"][1].replace(" + R", "\n+ R")
    new_data.loc[2, "Source"] = new_data["Source"][2].replace(" +", "\n+")
    new_data.loc[3, "Source"] = new_data["Source"][3].replace(" +", "\n+")
    new_data.loc[4, "Source"] = new_data["Source"][4].replace(" + S", "\n+ S")
    new_data.loc[6, "Source"] = new_data["Source"][6].replace(" Pl", "\nPl")

    labels: list[str] = list(new_data["Source"])
    labels_short: list[str] = [
        "Reg. BB",
        "SOA + SOIL\n+ RD",
        "Vehicles",
        "Const.+ Grills",
        "Ship, Trucks\nSS",
        "Agricultural",
        "TPP + Ind.",
    ]
    labels_short = [
        labels_short[j]
        + " "
        + new_data["PM2,5"].round(0).astype("int").astype("str")[j]
        + "%"
        for j in range(0, len(labels_short))
    ]
    values: np.ndarray[tuple[int], np.dtype[np.float64]
                       ] = new_data["PM2,5"].values
    fig: plt.Figure
    axs: plt.Axes

    cmap = mpl.colormaps["magma_r"]
    cmap = cmocean.tools.lighten(cc.colormaps["cet_glasbey"], 0.65)
    colors = cmap(np.linspace(0, 1, 256))
    fig, axs = plt.subplots(
        figsize=(8, 2.5), layout="constrained", ncols=2, width_ratios=[4, 1]
    )

    data2.plot.bar(
        stacked=True,
        color=colors[:],
        ax=axs[0],
        fontsize=7,
        legend=False,
        ylabel="Contribution (%)",
    )

    wedges, texts = axs[1].pie(values, startangle=45, colors=colors[:])

    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0)
    kw = dict(arrowprops=dict(arrowstyle="-"),
              bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1) / 2.0 + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = f"angle,angleA=0,angleB={ang}"
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        axs[1].annotate(
            labels_short[i],
            xy=(x, y),
            # xytext=(1.35*np.sign(x), 1.4*y),
            xytext=(1.35 * x, 1.4 * y),
            horizontalalignment=horizontalalignment,
            fontsize=7,
            ** kw,
        )
    handles, _ = axs[0].get_legend_handles_labels()
    print(labels)
    fig.legend(handles, labels, loc='lower center',
               ncols=7, fontsize=7, bbox_to_anchor=(0.5, -0.15))

    fig.savefig("PMF_contrib_stacked.png", dpi=300, bbox_inches="tight")
    plt.close()
    plt.show()


if __name__ == "__main__":
    main()
