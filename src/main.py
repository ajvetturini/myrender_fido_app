import os
import sys
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QUrl
from app import run_server
import atexit
from PyQt5.QtWebEngineWidgets import QWebEngineSettings
import multiprocessing
import time
import requests


# Global reference to the server process
server_process = None


# Create a multiprocessing target function for the server
def start_server():
    run_server()

def is_server_responsive():
    # Check if the server is responsive by sending  request
    try:
        response = requests.get("http://localhost:8080")
        return True
    except Exception:
        return False

def delete_files_in_directory(directory_path):
    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)
        if os.path.isfile(file_path):
            if '.DS_Store' not in file_path:
                os.remove(file_path)


# Function to be called when the application exits
def exit_handler():
    global server_process
    # This exit handler will kill the running server AND delete any files found within the flask_session or temp folders
    if server_process is not None and server_process.is_alive():
        # First delete files
        curDir = os.getcwd()
        dir1 = os.path.join(os.path.join(curDir, 'src'), 'flask_session')
        dir2 = os.path.join(os.path.join(curDir, 'src'), 'temp_stored_data')
        dir3 = os.path.join(os.path.join(curDir, 'src'), 'temp_uploaded_files')
        delete_files_in_directory(dir1)
        delete_files_in_directory(dir2)
        delete_files_in_directory(dir3)

        # Next end the server processes to close the application:
        server_process.terminate()
        server_process.join()


def main():
    global server_process
    # Start the server in a separate process within the "if __name__" block
    server_process = multiprocessing.Process(target=start_server)
    server_process.start()

    while not is_server_responsive():
        time.sleep(2)  # Sleep for 2 seconds and check again

    # Register the exit handler to clean up the server process
    atexit.register(exit_handler)

    # Create the PyQt5 app
    app = QApplication(sys.argv)

    class DashWindow(QMainWindow):
        def __init__(self, url):
            super(DashWindow, self).__init__()  # Call the parent class constructor correctly

            self.browser = QWebEngineView()
            self.browser.setUrl(QUrl(url))
            # Handling Downloads:
            self.browser.page().profile().downloadRequested.connect(self.on_downloadRequested)

            # Configure QWebEngineSettings
            web_settings = self.browser.settings()
            web_settings.WebAttribute.JavascriptEnabled = True
            web_settings.WebAttribute.JavascriptCanOpenWindows = True
            web_settings.WebAttribute.PluginsEnabled = True
            web_settings.WebAttribute.LocalStorageEnabled = True
            web_settings.WebAttribute.Antialiasing = True
            web_settings.WebAttribute.AutoLoadImages = True
            web_settings.WebAttribute.XSSAuditingEnabled = True
            web_settings.WebAttribute.ScrollAnimatorEnabled = True

            self.setCentralWidget(self.browser)
            self.setWindowTitle("AJ Project Name WIP")
            self.showMaximized()

        @QtCore.pyqtSlot("QWebEngineDownloadItem*")
        def on_downloadRequested(self, download):
            download_path = download.path()
            curSuffix = QtCore.QFileInfo(download_path).suffix()
            if curSuffix == 'json':
                old_path = 'MOSA_Input'
                suffix = 'json'
            else:
                old_path = 'MOSA_Output_Files'
                suffix = 'zip'

            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self, "Save File", old_path, "*." + suffix
            )
            if path:
                download.setPath(path)
                download.accept()


    dash_url = "http://localhost:8080"  # Replace with the actual URL of your Dash app
    dash_window = DashWindow(dash_url)
    dash_window.show()

    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
