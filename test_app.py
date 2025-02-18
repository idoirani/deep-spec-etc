
import dash
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = html.Div([
    dcc.Store(id='is-mobile-store', data=False),
    dcc.Interval(id='interval', interval=500, n_intervals=0, max_intervals=3),
    html.Div(id='output')
])

app.clientside_callback(
    """
    function(n_intervals) {
        console.log("n_intervals:", n_intervals);
        var isMobile = /Mobi|Android|iPhone|iPad/i.test(navigator.userAgent);
        console.log("Detected mobile:", isMobile);
        return isMobile;
    }
    """,
    Output('is-mobile-store', 'data'),
    Input('interval', 'n_intervals')
)

app.clientside_callback(
    """
    function(is_mobile) {
        return "Mobile flag is: " + is_mobile;
    }
    """,
    Output('output', 'children'),
    Input('is-mobile-store', 'data')
)

if __name__ == '__main__':
    app.run_server(debug=True)
