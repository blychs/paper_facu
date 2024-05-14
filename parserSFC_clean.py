import pandas as pd
import numpy as np


def readSFC(path2SFC: str) -> pd.DataFrame:
    """Parse *.SFC file"""
    names: list[str] = [
        "year",
        "month",
        "day",
        "julianDay",
        "hour",
        "H",
        "Ustar",
        "Wstar",
        "dtheta_dz",
        "PBLHc",
        "PBLHm",
        "L",
        "z0",
        "b0",
        "a0",
        "ws",
        "wd",
        "wind_h",
        "temp",
        "temp_h",
        "pp_code",
        "pp",
        "rh",
        "pres",
        "cloudCover",
        "es_adj",
        "subst",
    ]
    data: pd.DataFrame = pd.read_csv(path2SFC, skiprows=0, sep=r"\s+", names=names)
    data["year"] = data["year"] + 2000
    print(data["year"])
    data["date"] = pd.to_datetime(
        data[["year", "month", "day", "hour"]], format="%Y-%m-%d %H"
    )
    data["H"] = data["H"].replace(-999.0, np.nan)
    data["Ustar"] = data["Ustar"].replace(-9.0, np.nan)
    data["Wstar"] = data["Wstar"].replace(-9.0, np.nan)
    data["dtheta_dz"] = data["dtheta_dz"].replace(-9.0, np.nan)
    data["PBLHc"] = data["PBLHc"].replace(-999.0, np.nan)
    data["PBLHm"] = data["PBLHm"].replace(-999.0, np.nan)
    data["L"] = data["L"].replace(-99999.0, np.nan)
    data["z0"] = data["z0"].replace(-9.0, np.nan)
    data["b0"] = data["b0"].replace(-9.0, np.nan)
    data["a0"] = data["a0"].replace(-9.0, np.nan)
    data["ws"] = data["ws"].replace(999.0, np.nan)
    data["wd"] = data["wd"].replace(999.0, np.nan)
    data["temp"] = data["temp"].replace(-999.0, np.nan)
    data["pp"] = data["pp"].replace(-9.0, np.nan)
    data["rh"] = data["rh"].replace(999, np.nan)
    data["pres"] = data["pres"].replace(-999.0, np.nan)

    return data


def total_pblh(data: pd.DataFrame) -> pd.Series:
    """Calculate PBLH using L, from .SFC data"""

    data2: pd.DataFrame = data.copy()
    data2["maxPBLH"] = data2[["PBLHc", "PBLHm"]].max(axis=1)
    data2["PBLHLneg"] = data2["maxPBLH"].where(data2["L"] <= 0)
    data2["PBLHLpos"] = data2["PBLHm"].where(data2["L"] > 0)
    data2["PBLH"] = data2[["PBLHLneg", "PBLHLpos"]].sum(axis=1)
    data2["PBLH"] = data2["PBLH"].where(data2["L"].notna())

    return data2["PBLH"]


def ventilation_coef(data: pd.DataFrame) -> pd.Series:
    """Calculate ventilation coefficient using L, from existing data"""

    data2: pd.DataFrame = data.copy()
    data2["ventCoef"] = data2["PBLH"] * data2["ws"]

    return data2["ventCoef"]


def select_dates(data: pd.DataFrame, target_data: pd.DataFrame) -> pd.DataFrame:
    """Select data in pd.DataGrame data from dates present in target data"""

    data_out: pd.DataFrame = data.loc[
        data["date"].dt.date.isin(target_data["date"].dt.date)
    ]
    return data_out


def main() -> None:
    data_obs: pd.DataFrame = readSFC("OBSERVATORIO.SFC")
    dataEra5: pd.DataFrame = pd.read_csv("blh20192020CNEACent.txt")
    dataEra5["date"] = pd.to_datetime(dataEra5["date"])
    dataEra5 = dataEra5[["date", "blh"]]
    data_out: pd.DataFrame = pd.merge(data_obs, dataEra5, on="date")
    data_out = data_out.set_index(data_out["date"])
    data_out = data_out[["date", "ws", "temp", "rh", "pres", "blh"]]

    data_out.to_csv("test.csv")
    data_out_noon2noon = data_out.resample("24h", offset="12h").mean()

    target_data: pd.DataFrame = pd.read_excel("PMF_BA_fullv3.xlsx", sheet_name="CONC")
    data_out_noon2noon = select_dates(data_out_noon2noon, target_data)
    data_out_noon2noon["temp"] = data_out_noon2noon["temp"] - 273.15
    data_out_noon2noon[["ws", "temp", "rh", "pres", "blh"]].to_csv(
        "datos_meteo_blhera5.csv", float_format="%.2f"
    )


if __name__ == "__main__":
    main()
