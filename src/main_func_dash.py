import os

import dash
from dash import html

def create_dash_application(flask_app):
    #assets_folder = os.path.join(os.getcwd(), 'static')
    assets_folder = os.path.join(os.path.join(os.getcwd(), 'src'), 'static')
    dash_app = dash.Dash(
        server=flask_app,
        name='FIDO Results Dashboard',
        url_base_pathname='/',  # Used to be /results/
        assets_folder=assets_folder,
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
