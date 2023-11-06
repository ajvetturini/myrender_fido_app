from flask import Blueprint, render_template, flash, redirect, url_for, session, current_app, request
from werkzeug.utils import secure_filename
from flask_wtf import FlaskForm
from flask_wtf.file import FileField
from wtforms import SubmitField
from wtforms.validators import InputRequired
import pickle
import uuid
import os
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from .MultiObjective_Shape_Annealing import ShapeAnnealer

main = Blueprint('main', __name__)
ALLOWED_EXTENSIONS = {'aj1'}


class UploadForm(FlaskForm):
    file = FileField('Upload .aj1 File', validators=[InputRequired()])
    submit = SubmitField('Upload')

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def verify_uploaded_file(file):
    try:
        loaded_data = pickle.load(file.read())
        print(loaded_data)
        if isinstance(loaded_data, ShapeAnnealer):
            return True
        return False
    except Exception as e:
        print('ERROR: Invalid File.')
        return True

#@main.route('/', methods=['GET', 'POST'])
@main.route('/upload_results_file', methods=['GET', 'POST'])
def main_index():
    file_threshold = 100  # Max number of files we can store
    form = UploadForm()

    if form.validate_on_submit():
        file = form.file.data
        if file and allowed_file(file.filename):
            if verify_uploaded_file(file.stream):
                file.seek(0)  # Go to start
                # Save the uploaded file to my UPLOAD_FOLDER:
                unique_id_to_serialize = str(uuid.uuid4())
                filename = f"{unique_id_to_serialize}_{secure_filename(file.filename)}"
                upload_folder = os.path.join(current_app.config['UPLOAD_FOLDER'])
                number_of_files = len([entry for entry in os.listdir(upload_folder) if os.path.isfile(os.path.join(upload_folder, entry))])
                if number_of_files > file_threshold:
                    flash("Currently there are too many files stored on our server. These are cleaned out every hour, so please check back in an hour", "error")
                else:
                    # If there is enough space, we can go to the results dashboard!
                    file.save(os.path.join(current_app.config['UPLOAD_FOLDER'], filename))
                    session['user_id'] = unique_id_to_serialize
                    session['data_filepath'] = filename
                    if hasattr(file, 'close'):
                        file.close()  # Prevent memory leak...
                    return redirect('/results/')
            else:
                flash("Corrupted data (or invalid data object), please upload a .aj1 file.", "error")
        else:
            flash("Invalid file type. Please upload a .aj1 file.", "error")

    return render_template('home.html', active_page='homepage', form=form)

# Results dashboard redirection. OLD CODE, leaving here just for later just in case...?
"""@main.route('/results')
def results():
    #return redirect(url_for('main.results'))
    return render_template('github_n_links.html', active_page='other')"""


## ABOUT PAGE ROUTING:
@main.route("/about")
def about():
    return render_template('about.html', active_page='about')

## GITHUB PAGE ROUTING:
@main.route('/')
@main.route('/home')
def see_links():
    return render_template('github_n_links.html', active_page='other')

## INPUT CREATION PAGE ROUTING:
@main.route("/create_input_file")
def create_input_file():
    # If we are not rendering the simulation result, then we render the input form:
    ## This dictionary contains a list of all the min and max values to use:
    min_max_data = {
        'min_epochs': 10,
        'max_epochs': 500,
        'min_NT1': 500,
        'max_NT1': 5000,
        'min_NT2': 250,
        'max_NT2': 2500,
        'min_NBi': 500,
        'max_NBi': 5000,
        'min_Na': 100,
        'max_Na': 1500,
        'min_rb': 0.75,
        'max_rb': 0.99,
        'min_Nbi_LowerLim': 10,
        'max_Nbi_LowerLim': 50,
        'min_setsize': 10,
        'max_setsize': 25,
        'min_ri': 0.75,
        'max_ri': 0.99,
        'min_time_allowed': 5,
        'max_time_allowed': 180,
        'min_allowableT': 1e-12,
        'max_allowableT': 1e-5,
        'min_maxX': 40,
        'max_maxX': 100,
        'min_maxY': 40,
        'max_maxY': 100,
        'min_maxZ': 40,
        'max_maxZ': 100,
        'min_bps': 3000,
        'max_bps': 100000,
        'min_minEL': 39,
        'max_minEL': 105,
        'min_minAng': 8,
        'max_minAng': 20,
        'min_extend_rule_dist': 1,
        'max_extend_rule_dist': 5,
        'min_multiplier': 7,
        'max_multiplier': 15,
        'min_repulsion_const': 0.1,
        'max_repulsion_const': 100,
        'min_repulsion_dist': 0.5,
        'max_repulsion_dist': 20,
        'min_edges_at_node': 5,
        'max_edges_at_node': 10,
        'min_decimals': 5,
        'max_decimals': 10,
        'min_cooling_rate': 0.5,
        'max_cooling_rate': 0.99,
        'min_triki_T': 0.001,
        'max_triki_T': 100

    }
    ## Plotting "Dummy" Data to render a table:
    fig = make_subplots(rows=1, cols=1)
    fig.update_layout(title='Plotly Example')
    graphJSON = fig.to_json()

    return render_template('create_input_file.html', active_page='create_input_file', data=min_max_data,
                           graphJSON=graphJSON)

## VARIABLE DEFINITIONS:
@main.route("/parameter_definitions")
def parameter_definitions():
    return render_template('parameter_definitions.html', active_page='parameter_definitions')
