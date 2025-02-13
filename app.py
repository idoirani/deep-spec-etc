# -*- coding: utf-8 -*-


import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import parameters  # your parameters module
import ETC         # our refactored ETC module
import numpy as np 
# Initialize the Dash app.
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Define the layout.
app.layout = dbc.Container([
    # Title Frame
    dbc.Row(
        dbc.Col(
            dbc.Card([
                dbc.CardBody(
                    html.H1("MAST/DeepSpec SNR Calculator", className="text-center")
                )
            ]),
            width=6
        ),
        className="mb-4"
    ),
dbc.Row(
    dbc.Col(
        dbc.Card([
            dbc.CardBody([
                html.H5(
                    "This is our BETA version, and is still under development and testing (updated 13/02/2024)",
                    className="text-left",
                    style={"color": "red"}
                ),
                html.P(
                    "The calculations in this ETC are based on the as-built wavelength dependent resolution, and estimated end-to-end throughput of DeepSpec coupled to the MAST array, and assuming no slit losses (see https://ui.adsabs.harvard.edu/abs/2024SPIE13094E..53I/abstract). For comments and questions, please contact idoirani@gmail.com. ",
                    style={"color": "black", "textAlign": "left"}
                )
            ])
        ]),
        width=6
    ),
    className="mb-4"
),
    # Graph Frame
    dbc.Row(
        dbc.Col(
            dbc.Card([
                #dbc.CardHeader("Graph"),
                dbc.CardBody(
                    dcc.Graph(id="graph-output")
                )
            ]),
            width=6
        ),
        className="mb-4"
    ),
    # Parameter Card row
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
                                    {"label": "Spectra per hour", "value": "spec_per_hour"}
                                ],
                                value=parameters.calc_type
                            )
                        ], md=4),
                        dbc.Col([
                            html.Label("Sky Brightness (relative to dark)"),
                            dcc.Dropdown(
                                id="sky-brightness",
                                options=[
                                    {"label": "Dark", "value": 0},
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
                                    {"label": "[2,2]", "value": "[2,2]"}
                                ],
                                value='[1,1]'
                            )
                        ], md=4)
                    ], className="mb-3"),
                    
                    # SNR-specific parameters row
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
                                            {"label": "Pickels' Stellar Spectra", "value": "stellar"}
                                        ],
                                        value=parameters.spec_type
                                    )
                                ], md=4),
                                # Stellar spectrum options (conditionally visible via callback)
                                dbc.Col(
                                    html.Div(
                                        id="stellar-spectrum-options",
                                        children=[
                                            html.Label("Choose Stellar Type"),
                                            dcc.Dropdown(
                                                id="stellar-spectrum-dropdown",
                                                options=[{'label': key, 'value': key} for key in parameters.path_dic.keys()],
                                                value="G"
                                            )
                                        ],
                                        style={"display": "none"}  # hidden by default

                                    ),
                                    md=4
                                ),
                                # blackbody temperature (conditionally visible via callback)

                                dbc.Col(
                                    html.Div(
                                        id="bb-temp-container",
                                        children=[
                                            html.Label("Blackbody Temperature"),
                                            dcc.Input(id="bb-temp", type="number", value=parameters.bb_temp)
                                        ]
                                    ),
                                    md=4
                                ),
                            ], className="mb-2"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Exposure Time"),
                                    dcc.Input(id="t_snr", type="number", value=parameters.T_norm)
                                ], md=4), 

                                dbc.Col([
                                    html.Label("V Band Magnitude"),
                                    dcc.Input(id="ab-mag-renorm", type="number", value=parameters.AB_mag_renorm)
                                ], md=4),
                                dbc.Col([
                                    html.Label("Number of Telescopes in Group"),
                                    dcc.Input(id="n-tel-group", type="number", value=parameters.n_tel_group)
                                ], md=4)
                            ], className="mb-2"),

                        ]),
                        id = "snr-options-row",
                        className="mb-3"
                    ),
                    
                    # Limiting Magnitude-specific parameters row
                    dbc.Row(
                        dbc.Col([
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Number of Telescopes"),
                                    dcc.Input(id="n-tel-arr", type="text", value="1,4,10,20")
                                ], md=4),
                                dbc.Col([
                                    html.Label("Maximum Exposure Time"),
                                    dcc.Input(id="t_max", type="number", value=parameters.T_exp)
                                ], md=4),
                                dbc.Col([
                                    html.Label("Wavelength (Ang)"),
                                    dcc.Input(id="wl", type="number", value=parameters.wl)
                                ], md=4)
                            ], className="mb-2"),
                            dbc.Row([
                                dbc.Col([
                                    html.Label("SNR Limit"),
                                    dcc.Input(id="sigma_limit", type="number", value=parameters.sigma_limit)
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
                    
                    # Spectra per Hour-specific parameters row
                    dbc.Row(
                        dbc.Col([
                            dbc.Row([
                                dbc.Col([
                                    html.Label("Number of Telescopes in Full Array"),
                                    dcc.Input(id="n-tel", type="number", value=parameters.n_tel)
                                ], md=4),
                                dbc.Col([
                                    html.Label("Overhead (sec)"),
                                    dcc.Input(id="overhead", type="number", value=parameters.overhead_sec)
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
            width=6
        )
    )
], fluid=True)

# --- Callbacks ---

# Toggle visibility of Blackbody Temperature (only when spec-type is "bb")
@app.callback(
    Output("bb-temp-container", "style"),
    Input("spec-type", "value")
)
def toggle_bb_temp_options(spec_type):
    if spec_type == "bb":
        return {"display": "block", "margin-bottom": "1rem"}
    else:
        return {"display": "none"}

# Toggle visibility of Stellar Spectrum Options (only when spec-type is "stellar")
@app.callback(
    Output("stellar-spectrum-options", "style"),
    Input("spec-type", "value")
)
def toggle_stellar_options(spec_type):
    if spec_type == "stellar":
        return {"display": "block", "margin-bottom": "1rem"}
    else:
        return {"display": "none"}

# Toggle visibility of entire option rows based on calculation type.
@app.callback(
    Output("snr-options-row", "style"),
    Output("limmag-options-row", "style"),
    Output("spec-options-row", "style"),
    Input("calc-type", "value")
)
def toggle_calc_type_rows(calc_type):
    hidden = {"display": "none"}
    visible = {"display": "block"}
    if calc_type == "SNR":
        return visible, hidden, hidden
    elif calc_type == "limmag":
        return hidden, visible, hidden
    elif calc_type == "spec_per_hour":
        return hidden, hidden, visible
    else:
        return hidden, hidden, hidden

# --- Callback to Update the Plot ---
@app.callback(
    Output("graph-output", "figure"),
    Input("refresh-button", "n_clicks"),
    State("calc-type", "value"),
    State("spec-type", "value"),
    State("stellar-spectrum-dropdown", "value"),
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
def update_plot(n_clicks, calc_type, spec_type, stellar_type, bb_temp, ab_mag_renorm,
                n_tel, wl, n_tel_arr, n_tel_group, T_exp, sigma_limit, Type,
                delta_sky, binning, T_norm, overhead_sec, SNR_arr):
    import importlib
    importlib.reload(parameters)
    parameters.calc_type = calc_type
    parameters.spec_type = spec_type
    parameters.bb_temp = bb_temp
    parameters.AB_mag_renorm = ab_mag_renorm
    parameters.n_tel = n_tel
    parameters.wl = wl
    parameters.n_tel_arr = [int(x.strip()) for x in n_tel_arr.split(',')]
    parameters.n_tel_group = n_tel_group
    parameters.T_exp = T_exp
    parameters.sigma_limit = sigma_limit
    parameters.Type = Type
    parameters.delta_sky = delta_sky
    parameters.T_norm = T_norm
    parameters.overhead_sec = overhead_sec
    parameters.SNR_arr = [int(x.strip()) for x in SNR_arr.split(',')]
    if spec_type == "stellar":
        parameters.stellar_type = stellar_type
    else:
        parameters.stellar_type = None
    fig = ETC.run_ETC()
    return fig

if __name__ == "__main__":
    app.run_server(debug=True)