from shiny import App, render, ui

# Import modules for plot rendering
# %%
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
from scipy.stats import linregress, spearmanr, zscore
import statsmodels.api as sm
from funciones_pmfBA import mass_reconstruction, mass_reconstruction_mod, percentage_with_err
from funciones_pmfBA import estimation_om_oc, calculate_seasonal
from load_data import load_data


matrix, unc, meteo, gases, events = load_data('PMF_BA_full.xlsx', 'PMF_BA_full.xlsx',
                                               'gases_mean.csv', 'datos_meteo_obs_mean.csv',
                                               'BA_events.xlsx')


choices = list(matrix.keys())

app_ui = ui.page_fluid(
    ui.input_selectize("x", "Select x (single)", choices),
    ui.input_select("y", "Select y (single)", choices),
    ui.panel_main(
        ui.output_plot("scatter", width='600px', height='600px'),
    ),
)


def server(input, output, session):
    @output
    @render.plot(alt="A scatter plot")
    def scatter():
        fig, ax = plt.subplots(figsize=(20,20))
        ax.plot(matrix[input.x()], matrix[input.y()], 'o')
        ax.set_xlabel(input.x())
        ax.set_ylabel(input.y())
        print(ax.get_ylim()[0])

        



app = App(app_ui, server, debug=True)
