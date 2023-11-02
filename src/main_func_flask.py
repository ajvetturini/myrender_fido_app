"""
This will support the development of DNA Wireframe Modules using a simple and minimal web based application that can
then be hosted.
"""
from flask import Flask, current_app
from flask_session import Session
from datetime import timedelta
from views import main
import os


def create_app():
    myapp = Flask(__name__, template_folder='/Users/kodak/FlaskApp_FIDO_Host/src/templates', static_folder='/Users/kodak/FlaskApp_FIDO_Host/src/static')
    basedir = os.path.abspath(os.path.dirname(__file__))
    UPLOAD_FOLDER = os.path.join(basedir, 'temp_uploaded_files')
    myapp.config['SESSION_TYPE'] = 'filesystem'
    myapp.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
    myapp.config['SESSION_COOKIE_NAME'] = 'your_app_session'  # Set a session cookie name
    myapp.config['PERMANENT_SESSION_LIFETIME'] = timedelta(minutes=60)  # Set session expiration time
    myapp.secret_key = os.urandom(24)
    myapp.register_blueprint(main)
    Session(myapp)  # Store into the session

    return myapp


app = create_app()
