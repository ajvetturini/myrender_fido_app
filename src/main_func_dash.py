import dash
from dash import html

def create_dash_application(flask_app):
    dash_app = dash.Dash(
        server=flask_app,
        name='FIDO Results Dashboard',
        url_base_pathname='/',  # Used to be /results/
        assets_folder='/Users/kodak/FlaskApp_FIDO_Host/src/static',
        suppress_callback_exceptions=True
    )
    dash_app.layout = html.Div(
        children=[
            html.H1(children="FIDO Results Analysis"),
            html.Div(
                children=""
                         "Dashboard created via Dash"
            ),

        ]
    )
    return dash_app
