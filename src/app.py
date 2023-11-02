from main_func_flask import app as flask_app
from main_func_dash import create_dash_application  # Import your Dash app from dash_app.py
from dash_results_board import update_dashboard
import threading
import time
import os
from flaskwebgui import FlaskUI

def cleanup_task():
    while True:
        temp_folder = os.path.join(os.getcwd(), 'temp_uploaded_files')
        if not os.path.exists(temp_folder):
            os.makedirs(temp_folder)
        current_time = time.time()
        max_age = 3600  # This is in SECONDS, and the max age for a file is 1 hour
        # This function will loop through all stored data and uploaded files and delete them.
        for file_name in os.listdir(temp_folder):
            file_path = os.path.join(temp_folder, file_name)
            if os.path.isfile(file_path) and (current_time - os.path.getctime(file_path) > max_age):
                # Remove the file
                os.remove(file_path)
        # After this, we go to sleep for the set amount of age time:
        time.sleep(max_age)


app = create_dash_application(flask_app=flask_app)
server = app.server
# Call the update_dashboard function to set up the initial layout
update_dashboard(app)
# Create a simple cleanup watcher:
#cleanup_thread = threading.Thread(target=cleanup_task)
#cleanup_thread.daemon = True
#cleanup_thread.start()

def run_server():
    app.run_server(host='0.0.0.0', port=8080, debug=False)


if __name__ == '__main__':
    app.run_server(host='0.0.0.0', port=8080, debug=False)

