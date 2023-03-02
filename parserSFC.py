import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def readSFC(path2SFC):
    names = ["year", "month", "day", "julianDay", "hour",
            "H", "Ustar", "Wstar", "dtheta_dz", "PBLHc", "PBLHm", "L",
            "z0", "b0", "a0", "ws", "wd", "wind_h", "temp", "temp_h",
            "pp_code", "pp", "rh", "pres", "cloudCover", "es_adj",
            "subst"]
    data = pd.read_csv(path2SFC, skiprows=1, delim_whitespace=True,
                       names=names)
    data["year"] = data["year"] + 2000
    data["date"] = pd.to_datetime(data[["year", "month", "day", "hour"]])
    data = data.set_index(data["date"])
    data["H"].replace(-999.0, np.nan, inplace=True)
    data["Ustar"].replace(-9.0, np.nan, inplace=True)
    data["Wstar"].replace(-9.0, np.nan, inplace=True)
    data["dtheta_dz"].replace(-9.0, np.nan, inplace=True)
    data["PBLHc"].replace(-999.0, np.nan, inplace=True)
    data["PBLHm"].replace(-999.0, np.nan, inplace=True)
    data["L"].replace(-99999.0, np.nan, inplace=True)
    data["z0"].replace(-9.0, np.nan, inplace=True)
    data["b0"].replace(-9.0, np.nan, inplace=True)
    data["a0"].replace(-9.0, np.nan, inplace=True)
    data["ws"].replace(999.0, np.nan, inplace=True)
    data["wd"].replace(999.0, np.nan, inplace=True)
    data["temp"].replace(-999.0, np.nan, inplace=True)
    data["pp"].replace(-9.0, np.nan, inplace=True)
    data["rh"].replace(999, np.nan, inplace=True)
    data["pres"].replace(-999.0, np.nan, inplace=True)

    return data

def ventilation_coef(data):
    """Calculate ventilation coefficient using
    L"""
    data2 = data.copy()
    data2["maxPBLH"] = data2[["PBLHc", "PBLHm"]].max(axis=1)
    data2["VCLneg"] = data2[["ws", "maxPBLH"]].product(axis=1).where(
                                               data2["L"] <= 0)
    data2["VCLpos"] = data2[["ws", "PBLHm"]].product(axis=1).where(
                                             data2["L"] > 0)
    data2["ventCoef"] = data2[["VCLneg", "VCLpos"]].sum(axis=1)
    data2["ventCoef"] = data2["ventCoef"].where(data2["L"].notna())
    
    return data2["ventCoef"]

def diurnal_cycle(data):
    """ Calculate the diurnal cycle of hourly data with H as hour
    column"""

    data_diurnal = data.groupby('hour').mean()
    data_diurnal_err = data.groupby('hour').std()
    return data_diurnal, data_diurnal_err

def main():
    dataObs = readSFC("UAESTIMATOR/OBSERVATORIO.SFC")
    dataObs["ventCoef"] = ventilation_coef(dataObs)
    dataObs_noon2noon = dataObs.resample('24H', offset='12H').mean()
    dataObs_noon2noon.to_csv('dataObs_noon2noon.csv')
    # dataObsMed = readSFC("ORTUZAR/USANDO RAOB EZEIZA/OBSERVATORIO.SFC")
    # dataObsMed["ventCoef"] = ventilation_coef(dataObsMed)
    dataEzeEst = readSFC("UAESTIMATOR/EZEIZA.SFC")
    dataEzeEst["ventCoef"] = ventilation_coef(dataEzeEst)
    dataEzeEst_diurnal, dataEzeEst_diurnal_err = diurnal_cycle(dataEzeEst)
    dataEzeAer = readSFC("AERMET/EZEIZA.SFC")
    dataEzeAer["ventCoef"] = ventilation_coef(dataEzeAer)

    dataObs_diurnal, dataObs_diurnal_err = diurnal_cycle(dataObs)


 #   dataEzeEst = dataEzeEst[dataEzeEst.index.isin(dataObs.index)]
    print(dataObs_diurnal[["L", "PBLHc", "PBLHm", "ws", "ventCoef"]])
    print(dataObs_diurnal_err[["L", "PBLHc", "PBLHm", "ws", "ventCoef"]])


    plt.style.use('seaborn-v0_8-paper')
#    plt.style.use('ggplot')

#    fig, ax = plt.subplots(figsize=(10,10))
#    ax.scatter(dataEze["PBLHm"], dataObs["PBLHm"])
#    ax.set_yscale("log")
#    ax.set_xscale("log")
#    ax.set_xlabel("PBLHm Ezeiza (m)")
#    ax.set_ylabel("PBLHm Observatorio (m)")
#    plt.show()
#    
#    fig, ax = plt.subplots(figsize=(10,10))
#    ax.scatter(dataEze["PBLHc"], dataObs["PBLHc"])
#    ax.plot([0,10e3], [0,10e3], color='k')
#    ax.set_yscale("log")
#    ax.set_xscale("log")
#    ax.set_xlabel("PBLHc Ezeiza (m)")
#    ax.set_ylabel("PBLHc Observatorio (m)")
#    plt.show()
#
#    fig, ax = plt.subplots(figsize=(20,10))
#    ax.plot(dataEze.index[500:596], dataEze['PBLHm'][500:596], '.', label='PBLHm')
#    ax.plot(dataEze.index[500:596], dataEze['PBLHc'][500:596], '.', label='PBLHc')
#    ax.set_xlabel('Datetime')
#    ax.set_ylabel('PBLH height (m)')
#    ax.legend()
#    plt.show()

    fig, ax = plt.subplots(figsize=(5,5))
    ax.errorbar(dataEzeEst_diurnal.index, dataEzeEst_diurnal['PBLHc'],
                yerr=dataEzeEst_diurnal_err['PBLHc'], capsize=3, capthick=2, label='PBLH c', zorder=1)
    ax.errorbar(dataEzeEst_diurnal.index, dataEzeEst_diurnal['PBLHm'],
                yerr=dataEzeEst_diurnal_err['PBLHm'], capsize=3, capthick=2, label='PBLH m', zorder=2)
    ax.plot(dataEzeEst_diurnal.index, dataEzeEst_diurnal['PBLHm'].where(dataEzeEst_diurnal["L"] < 0), 'Xk',
            label='L $\leq$ 0', zorder=3)
    ax.plot(dataEzeEst_diurnal.index, dataEzeEst_diurnal['PBLHm'].where(dataEzeEst_diurnal["L"] > 0), 'Xr',
            label='L $>$ 0', zorder=4)
    ax.axhline(0)
    ax.legend()
    fig.savefig('diurnal_pblh_ezeizaEstimator.png')
    plt.close()

    fig, ax = plt.subplots(figsize=(7,7))
    ax.scatter(dataEzeEst[dataEzeEst.index.isin(dataObs.index)]["PBLHc"], dataObs["PBLHc"])
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("PBLHc Ezeiza UA Est (m)", size=20)
    ax.set_ylabel("PBLHc Obs UA Est (m)", size=20)
    fig.savefig('PBLHc_EzeUAEst_ObsUAEst.png')
    plt.close()
    

    fig, ax = plt.subplots(figsize=(7,7))
    ax.scatter(dataEzeEst[dataEzeEst.index.isin(dataObs.index)]["PBLHm"], dataObs["PBLHm"])
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("PBLHm Eze UA Est (m)", size=20)
    ax.set_ylabel("PBLHm Obs UA Est (m)", size=20)
    fig.savefig('PBLHm_EzeUAest_ObsUAEst.png')
    plt.close()

    fig, ax = plt.subplots(figsize=(7,7))
    ax.scatter(dataEzeEst[dataEzeEst.index.isin(dataEzeAer.index)]["PBLHm"],
               dataEzeAer[dataEzeAer.index.isin(dataEzeEst.index)]["PBLHm"])
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("PBLHm Eze UA Est (m)", size=20)
    ax.set_ylabel("PBLHm Eze Aermet (m)", size=20)
    fig.savefig('PBLHm_Eze_Eze_EstAer.png')
    plt.close()

    fig, ax = plt.subplots(figsize=(7,7))
    ax.scatter(dataEzeEst[dataEzeEst.index.isin(dataEzeAer.index)]["PBLHc"],
               dataEzeAer[dataEzeAer.index.isin(dataEzeEst.index)]["PBLHc"])
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("PBLHc Eze UA Est (m)", size=20)
    ax.set_ylabel("PBLHc Eze Aermet (m)", size=20)
    fig.savefig('PBLHc_Eze_Eze_EstAer.png')
    plt.close()

    plt.show()


if __name__ == '__main__':
    main()