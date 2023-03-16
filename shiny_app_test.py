from shiny import App, render, ui, reactive

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
    ui.input_selectize("y", "Select y (single)", choices),
    ui.layout_sidebar(
        ui.panel_sidebar(
                ui.tags.p("Change valid xmax and ymax"),
                ui.input_slider("xlim", "Xlim", min=0, max=75, value=10),
                ui.input_slider("ylim", "Ylim", min=0, max=75, value=10)
        ),
    ui.panel_main(
        ui.output_plot("scatter", width='600px', height='600px'),
    ),
    )
)


def server(input, output, session):
    @output
    @render.plot(alt="A scatter plot")
    def scatter():
        x = matrix[input.x()]
        x = x.where(x <= input.xlim())
        y = matrix[input.y()]
        y = y.where(y <= input.ylim())

        mask = ~np.isnan(x) & ~np.isnan(y)
        slope, intercept, r, p, se = linregress(x[mask], y[mask])
        fig, ax = plt.subplots(figsize=(20,20))
        ax.plot(x, y, 'o')
        ax.plot(x, intercept + slope * x,
                 label=f'{slope:.2g} [{input.x()}] + {intercept:.2g} \nR$^2$ = {r**2:.2g}\np-val = {p:.2g}')
        ax.set_xlabel(input.x())
        ax.set_ylabel(input.y())
    #    ax.set_xlim(right=input.xmax())
        ax.legend()

        print(ax.get_ylim()[0])
    
    @reactive.Effect
    def _():
        #val = np.ceil(matrix[input.x()].max())
        # Check for the max value and increase by 1% to ensure no floatin point error
        valx = matrix[input.x()].max() * 1.01 
        valy = matrix[input.y()].max() * 1.01 
        ui.update_slider(
            "xlim", value=valx, min=0, max=valx, step=valx/15
        )
        ui.update_slider(
            "ylim", value=valy, min=0, max=valy, step=valy/15  
        )

        



app = App(app_ui, server, debug=True)
