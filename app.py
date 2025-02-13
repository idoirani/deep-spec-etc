# -*- coding: utf-8 -*- 

"""
This Dash application displays an Exposure Time Calculator (ETC) for the MAST/DeepSpec instrument.
It consists of several "frames" (cards) that contain the title, an informational notice, the plot,
and parameter controls. The parameter inputs are grouped into sections that are conditionally
displayed based on the selected calculation type (SNR, Limiting Magnitude, or Spectra per Hour).

The application is responsive and includes mobile-friendly meta tags.
"""

import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import parameters  # Module containing default parameters
import ETC         # Module that runs the ETC calculations and returns a Plotly figure
import numpy as np 
# Initialize the Dash app with Bootstrap CSS and mobile-friendly meta tags.
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP],
                meta_tags=[{'name': 'viewport', 'content': 'width=device-width, initial-scale=1'}])

# Define the layout.
app.layout = dbc.Container([
    # Title Frame: Displays the main title in a card.
    dbc.Row(
        dbc.Col(
            dbc.Card([
                dbc.CardBody(
                    html.H1("MAST/DeepSpec SNR Calculator", className="text-center")
                )
            ]),
            width=8
        ),
        className="mb-4"
    ),
    # Informational Frame: Displays a smaller red notice and an explanation in black.

    dbc.Row(
        dbc.Col(
            dbc.Card([
                dbc.CardBody([
                    html.H5(
                        "This is our BETA version, and is still under development and testing (updated 13/02/2024)",
                        className="text-left",
                        style={"color": "red"}
                    ),
                    # Explanation text in black on a new line.

                    html.P(
                        "The calculations in this ETC are based on the as-built wavelength dependent resolution, and estimated end-to-end throughput of DeepSpec coupled to the MAST array, and assuming no slit losses (see https://ui.adsabs.harvard.edu/abs/2024SPIE13094E..53I/abstract). For comments and questions, please contact idoirani@gmail.com. ",
                        style={"color": "black", "textAlign": "left"}
                    )
                ])
            ]),
            width=8
        ),
        className="mb-4"
    ),
    # Graph Frame: Displays the Plotly graph in its own card.
    dbc.Row(
        dbc.Col(
            dbc.Card([
                dbc.CardBody(
                    dcc.Graph(id="graph-output", config={'responsive': True})
                )
            ]),
            width=8
        ),
        className="mb-4"
    ),
    # Parameter Card Frame: Contains all parameter inputs grouped in several sections.
    dbc.Row(
        dbc.Col(
            dbc.Card([
                dbc.CardHeader("Observing Conditions and Source Parameters"),
                dbc.CardBody([
                    # Shared row (always visible)
                    dbc.Row([
                        dbc.Col([
                            html.Label("Calculation Type"),
                            dcc.Dropdown(
                                id="calc-type",
                                options=[
                                    {"label": "SNR", "value": "SNR"},
                                    {"label": "Limiting magnitude", "value": "limmag"},
                                    {"label": "Spectra per hour (this might take a moment)", "value": "spec_per_hour"}
                                ],
                                value=parameters.calc_type
                            )
                        ], md=4),
                        dbc.Col([
                            html.Label("Sky Brightness (relative to dark)"),
                            dcc.Dropdown(
                                id="sky-brightness",
                                options=[
                                    {"label": "Dark (20.9 mag/arcsec^2)", "value": 0},
                                    {"label": "Dark+2", "value": -2},
                                    {"label": "Dark+4", "value": -4},
                                    {"label": "Dark+6", "value": -6}
                                ],
                                value=-2
                            )
                        ], md=4),
                        dbc.Col([
                            html.Label("Binning"),
                            dcc.Dropdown(
                                id="binning",
                                options=[
                                    {"label": "[1,1]", "value": "[1,1]"},
                                    {"label": "[2,1]", "value": "[2,1]"},
                                    {"label": "[3,1]", "value": "[3,1]"},
                                    {"label": "[2,2]", "value": "[2,2]"},
                                    {"label": "[2,3]", "value": "[2,3]"},
                                    {"label": "[3,3]", "value": "[3,3]"},
                          ],
                                value='[1,1]'
                            )
                        ], md=4)
                    ], className="mb-3"),
                    
                    # SNR-Specific Parameters Row (only visible when calc-type == "SNR")
                    dbc.Row(
                        dbc.Col([
                            # Spectrum Type dropdown
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Spectrum Type"),
                                    dcc.Dropdown(
                                        id="spec-type",
                                        options=[
                                            {"label": "Blackbody", "value": "bb"},
                                            {"label": "Flat in f_nu", "value": "flat"},
                                            {"label": "Pickels' Stellar Spectra", "value": "stellar"},
                                            {"label": "Koester WD Models", "value": "WD"}

                                        ],
                                        value=parameters.spec_type
                                    )
                                ], md=4),
                                # Stellar spectrum options: Shown only when spec-type == "stellar"
                                dbc.Col(
                                    html.Div(
                                        id="stellar-spectrum-options",
                                        children=[
                                            html.Label("Choose Stellar Type"),
                                            dcc.Dropdown(
                                                id="stellar-spectrum-dropdown",
                                                options=[{'label': key, 'value': key} for key in parameters.stellar_path_dic.keys()],
                                                value="G"
                                            )
                                        ],
                                        style={"display": "none"}  # hidden by default

                                    ),
                                    md=4
                                ),
                                # Blackbody Temperature: Shown only when spec-type == "bb"

                                dbc.Col(
                                    html.Div(
                                        id="bb-temp-container",
                                        children=[
                                            html.Label("Blackbody Temperature"),
                                            dcc.Input(id="bb-temp", type="number", value=parameters.bb_temp, min=1000,max = 1000000, step = 1)
                                        ]
                                    ),
                                    md=4
                                ),
                                # Teff options: Shown only when spec-type == "WD" (conditional callback not shown here)

                                dbc.Col(
                                    html.Div(
                                        id="Teff-options",
                                        children=[
                                            html.Label("Teff"),
                                            dcc.Dropdown(
                                                id="Teff-dropdown",
                                                options=[{'label': T, 'value': T} for T in np.unique(parameters.Teff)],
                                                value="G"
                                            )
                                        ],
                                        style={"display": "none"}  # hidden by default

                                    ),
                                    md=4
                                ),
                                # logg options: Shown only when spec-type == "WD" (conditional callback not shown here)

                                dbc.Col(
                                    html.Div(
                                        id="logg-options",
                                        children=[
                                            html.Label("logg"),
                                            dcc.Dropdown(
                                                id="logg-dropdown",
                                                options=[{'label': Lg, 'value': Lg} for Lg in np.unique(parameters.logg)],
                                                value="G"
                                            )
                                        ],
                                        style={"display": "none"}  # hidden by default

                                    ),
                                    md=4
                                ),

                            ], className="mb-2"),
                            # Second row within SNR options: Exposure Time, V Band Magnitude, and number of telescopes.

                            dbc.Row([
                                dbc.Col([
                                    html.Label("Exposure Time"),
                                    dcc.Input(id="t_snr", type="number", value=parameters.T_norm, min=1,max = 2000, step = 1)
                                ], md=4), 

                                dbc.Col([
                                    html.Label("V Band Magnitude"),
                                    dcc.Input(id="ab-mag-renorm", type="number", value=parameters.AB_mag_renorm, min=5,max = 23)
                                ], md=4),
                                dbc.Col([
                                    html.Label("Number of Telescopes in Group"),
                                    dcc.Input(id="n-tel-group", type="number", value=parameters.n_tel_group, min=1,max = 20, step = 1)
                                ], md=4)
                            ], className="mb-2"),

                        ]),
                        id = "snr-options-row",
                        className="mb-3"
                    ),
                    
                    # Limiting Magnitude-Specific Parameters Row (only visible when calc-type == "limmag")
                    dbc.Row(
                        dbc.Col([
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Number of Telescopes"),
                                    dcc.Input(id="n-tel-arr", type="text", value="1,4,10,20")
                                ], md=4),
                                dbc.Col([
                                    html.Label("Maximum Exposure Time"),
                                    dcc.Input(id="t_max", type="number", value=parameters.T_exp, min=100,max = 2000, step = 1)
                                ], md=4),
                                dbc.Col([
                                    html.Label("Wavelength (Ang)"),
                                    dcc.Input(id="wl", type="number", value=parameters.wl, min=3700,max = 9000, step = 1)
                                ], md=4)
                            ], className="mb-2"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("SNR Limit"),
                                    dcc.Input(id="sigma_limit", type="number", value=parameters.sigma_limit, min=5, max = 100, step = 1)
                                ], md=4),
                                dbc.Col([
                                    html.Label("SNR Calculation Type"),
                                    dcc.Dropdown(
                                        id="snr-type",
                                        options=[
                                            {"label": "per pixel", "value": "per pixel"},
                                            {"label": "per resolution element", "value": "per element"}
                                        ],
                                        value=parameters.Type
                                    )
                                ], md=4)
                            ])
                        ]),
                        id = "limmag-options-row",
                        className="mb-3"
                    ),
                    
                    # Spectra per Hour-Specific Parameters Row (only visible when calc-type == "spec_per_hour")
                    dbc.Row(
                        dbc.Col([
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Number of Telescopes in Full Array"),
                                    dcc.Input(id="n-tel", type="number", value=parameters.n_tel, min=1,max = 20, step = 1)
                                ], md=4),
                                dbc.Col([
                                    html.Label("Overhead (sec)"),
                                    dcc.Input(id="overhead", type="number", value=parameters.overhead_sec, min=1,max = 1000, step = 1)
                                ], md=4),
                                dbc.Col([
                                    html.Label("SNR Array"),
                                    dcc.Input(id="snr-arr", type="text", value="10,20,50")
                                ], md=4)
                            ])
                        ]),
                        id="spec-options-row",
                        className="mb-3"
                    ),
                    
                    dbc.Row(
                        dbc.Col(
                            html.Button("Refresh Plot", id="refresh-button", n_clicks=0),
                            width={"size": 4, "offset": 4}
                        )
                    )
                    
                ])
            ]),
            width=8
        )
    )
], fluid=True)

# --- Callbacks ---

# Callback: Toggle visibility of Blackbody Temperature container.
@app.callback(
    Output("bb-temp-container", "style"),
    Input("spec-type", "value")
)
def toggle_bb_temp_options(spec_type):
    # Show Blackbody Temperature input only when spec-type is "bb"
    if spec_type == "bb":
        return {"display": "block", "margin-bottom": "1rem"}
    else:
        return {"display": "none"}

# Callback: Toggle visibility of Stellar Spectrum Options container.
@app.callback(
    Output("stellar-spectrum-options", "style"),
    Input("spec-type", "value")
)
def toggle_stellar_options(spec_type):
    # Show Stellar Spectrum dropdown only when spec-type is "stellar"
    if spec_type == "stellar":
        return {"display": "block", "margin-bottom": "1rem"}
    else:
        return {"display": "none"}

# Toggle visibility of Teff Spectrum Options (only when spec-type is "WD")
@app.callback(
    Output("Teff-options", "style"),
    Input("spec-type", "value")
)
def toggle_Teff_options(spec_type):
    # Show Teff dropdown only when spec-type is "WD"
    if spec_type == "WD":
        return {"display": "block", "margin-bottom": "1rem"}
    else:
        return {"display": "none"}
# Toggle visibility of logg Spectrum Options (only when spec-type is "WD")
@app.callback(
    Output("logg-options", "style"),
    Input("spec-type", "value")
)
def toggle_logg_options(spec_type):
    # Show logg dropdown only when spec-type is "WD"
    if spec_type == "WD":
        return {"display": "block", "margin-bottom": "1rem"}
    else:
        return {"display": "none"}
# Callback: Toggle visibility of entire option rows based on Calculation Type.
@app.callback(
    Output("snr-options-row", "style"),
    Output("limmag-options-row", "style"),
    Output("spec-options-row", "style"),
    Input("calc-type", "value")
)
def toggle_calc_type_rows(calc_type):
    hidden = {"display": "none"}
    visible = {"display": "block"}
    # Display the appropriate row based on selected calc-type.
    if calc_type == "SNR":
        return visible, hidden, hidden
    elif calc_type == "limmag":
        return hidden, visible, hidden
    elif calc_type == "spec_per_hour":
        return hidden, hidden, visible
    else:
        return hidden, hidden, hidden

# Callback: Update the Plot when the Refresh Plot button is clicked.
@app.callback(
    Output("graph-output", "figure"),
    Input("refresh-button", "n_clicks"),
    State("calc-type", "value"),
    State("spec-type", "value"),
    State("stellar-spectrum-dropdown", "value"),
    State("Teff-dropdown", "value"),
    State("logg-dropdown", "value"),
    State("bb-temp", "value"),
    State("ab-mag-renorm", "value"),
    State("n-tel", "value"),
    State("wl", "value"),
    State("n-tel-arr", "value"),
    State("n-tel-group", "value"),
    State("t_max", "value"),
    State("sigma_limit", "value"),
    State("snr-type", "value"),
    State("sky-brightness", "value"),
    State("binning", "value"),
    State("t_snr", "value"),
    State("overhead", "value"),
    State("snr-arr", "value")
)
def update_plot(n_clicks, calc_type, spec_type, stellar_type,Teff, logg, bb_temp, ab_mag_renorm,
                n_tel, wl, n_tel_arr, n_tel_group, T_exp, sigma_limit, Type,
                delta_sky, binning, T_norm, overhead_sec, SNR_arr):
    """
    Updates the global parameters based on user inputs and runs the ETC calculation.
    Returns a Plotly figure generated by ETC.run_ETC().
    """
    import importlib
    importlib.reload(parameters)
    parameters.calc_type = calc_type
    parameters.spec_type = spec_type
    parameters.bb_temp = bb_temp
    parameters.AB_mag_renorm = ab_mag_renorm
    parameters.n_tel = n_tel
    parameters.wl = wl
    # Convert the comma-separated telescope numbers to a list of integers.
    parameters.n_tel_arr = [int(x.strip()) for x in n_tel_arr.split(',')]
    parameters.n_tel_group = n_tel_group
    parameters.T_exp = T_exp
    parameters.sigma_limit = sigma_limit
    parameters.Type = Type
    parameters.delta_sky = delta_sky
    parameters.Sky_brightness_surface_den = parameters.best_sky_brightness_surface_den + parameters.delta_sky 

    parameters.T_norm = T_norm
    parameters.overhead_sec = overhead_sec
    # Convert the comma-separated SNR values to a list of integers.
    parameters.SNR_arr = [int(x.strip()) for x in SNR_arr.split(',')]
    # Convert binning string "[x,y]" to a list of integers.
    parameters.binning = [int(x) for x in binning.strip('[]').split(',')]

    # Set the stellar type if applicable.
    if spec_type == "stellar":
        parameters.stellar_type = stellar_type
    else:
        parameters.stellar_type = None
    # For WD models, set Teff and logg (if spec_type is 'WD')

    if spec_type == 'WD':
        parameters.Teff = Teff
        parameters.logg = logg
    else:
        parameters.Teff = None
        parameters.logg = None
    # Run the ETC calculations and return the resulting figure.

    fig = ETC.run_ETC()
    return fig

if __name__ == "__main__":
    app.run_server(debug=True)
    #app.run_server(debug=True, host='0.0.0.0', port=8091)