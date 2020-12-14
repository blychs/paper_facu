# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.7.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
from pandas.api.types import is_numeric_dtype
import xarray as xr


def corr_matrix(df, minimum=-1, maximum=1):
    return df.corr().where(df.corr() >= minimum).where(df.corr() <= maximum)


def test_lognormality(data_array, treshold=0):
    return stats.normaltest(np.log(data_array.where(data_array > treshold)), axis=0, nan_policy='omit').pvalue


def graph_all_corr(df):
    for i in df:
        if is_numeric_dtype(df[i]):
            for j in df:
                if is_numeric_dtype(df[j]):
                    idx = np.isfinite(df[i]) & np.isfinite(df[j])
                    coef = np.polyfit(df[i][idx],df[j][idx],1)
                    print(str(coef[0]), 'x +', str(coef[1]))
                    poly1d_fn = np.poly1d(coef)
                    correlation_fig = plt.plot(df[i],df[j], 'yo', df[i], poly1d_fn(df[i]), 'k')
                    plt.legend(iter(correlation_fig), ('Values', 'Linear fit'))
                    #plt.plot(df[i], df[i], 'r--')
                    plt.xlabel(i)
                    plt.ylabel(j)
                    plt.savefig('corr_' + i + '_' + j + '.png')
                    plt.show()
