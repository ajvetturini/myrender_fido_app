from dash import Dash, dcc, html, callback, callback_context
from dash.dependencies import Input, Output, State
from MultiObjective_Shape_Annealing import ShapeAnnealer
from Ply_Graphs_3D import DesignSpace3D
from tacoxDNA.src.CanDo_oxDNA import convert_cndo_to_oxdna_in_python
from MOSA_wireframe_generator import find_min_and_maxes
import plotly.graph_objs as go
import shutil
import random
import string
import os
import subprocess
import platform
import json
import pickle
from flask import session, current_app


def return_read_data():
    filename = session['data_filepath']
    filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
    with open(filepath, 'rb') as f:
        data = pickle.load(f)
    return data

def write_out_cando_cadnano(filepath, plypath, edge_length, dna_seq):
    writepath = os.path.join(filepath, 'caDNAno_cando_oxDNA_files')
    os.makedirs(writepath)
    DAED_MASTER_PATH = os.path.join(os.getcwd(), 'DAEDALUS')
    M13_PATH = os.path.join(DAED_MASTER_PATH, 'M13.txt')

    # Read in DNA sequence and determine if we need to save a txt file of the input scaffold:
    my_seq = json.loads(dna_seq)
    if my_seq == "NONE":
        SEQ_PATH = M13_PATH
        shutil.copy(M13_PATH, filepath)
    else:
        # If the user specified a DNA Sequence then we must save out to a txt file to use within DAEDALUS.
        custom_scaffold_path = os.path.join(filepath, 'custom_sequence.txt')
        with open(custom_scaffold_path, 'w') as file:
            file.write(my_seq)
        SEQ_PATH = custom_scaffold_path

    if platform.system() == "Windows":
        daed_path = os.path.join(DAED_MASTER_PATH, 'DAEDALUS2.exe')
        command = f"{daed_path} {writepath} {plypath} {SEQ_PATH} 1 2 0 {edge_length} 0.0 s"
    else:
        daed_path = os.path.join(DAED_MASTER_PATH, 'DAEDALUS2')
        command = f"{daed_path} {writepath} {plypath} {SEQ_PATH} 1 2 0 {edge_length} 0.0 s"

    subprocess.run(command, shell=True, text=True, check=True, stdout=subprocess.PIPE)  # Runs DAEDALUS 2 algorithm
    # Clean up the files from DAEDALUS:
    CNDO_FILEPATH = None
    rename_mapping = {
        '.cndo': 'GeneratedDesign_cando.cndo',
        '.csv': 'GeneratedDesign_scaffold_and_staple_sequences.csv',
        '.json': 'GenereatedDesign_caDNAno.json'
    }
    num_files_created = len([entry for entry in os.listdir(writepath) if os.path.isfile(os.path.join(writepath, entry))])
    # We check for the value 2 because if DAEDALUS fails due to invalid scaffold length, it always creates 2 bild files.
    if num_files_created == 2:
        file_path = os.path.join(writepath, 'ERROR_README.txt')
        with open(file_path, 'w') as f:
            f.write("DAEDALUS2 found that the required scaffold was longer than your input sequence (or longer than M13) and therefore an output was NOT generated. \n")
            f.write("This may be due to the way MOSA approximated the scaffold length (as it should prevent generating designs with longer scaffold lengths than the input length) and will be fixed in future versions! \n")
            f.write("Also, please check the 'max_scaffold_length' parameter in your MOSA process as that could be the source of the issue (you may be accidentally allowing longer scaffolds to be used than intended) \n")
        # Remove the 2 bild files afterwords
        for filename in os.listdir(writepath):
            file_path = os.path.join(writepath, filename)
            # Check if the file has a ".bild" extension and delete it
            if filename.endswith(".bild"):
                os.remove(file_path)
    else:
        for filename in os.listdir(writepath):
            file_path = os.path.join(writepath, filename)
            # Check if the file has a ".bild" extension and delete it
            if filename.endswith(".bild"):
                os.remove(file_path)
            # Rename some files:
            elif filename.endswith((".cndo", ".csv", ".json")):
                ext = os.path.splitext(filename)[1]
                new_file_path = None  # Initialize
                if ext in rename_mapping:
                    new_filename = rename_mapping[ext]
                    new_file_path = os.path.join(writepath, new_filename)
                    os.rename(file_path, new_file_path)

                # Also track for cndo extension:
                if filename.endswith(".cndo"):
                    CNDO_FILEPATH = new_file_path

    # Copy and paste the sequence text file to our save directory:


    # After cleaning out the folder, we simply pass the cando file_path to tacoxdna to convert and create oxDNA
    if CNDO_FILEPATH is not None:
        convert_cndo_to_oxdna_in_python(cndo_filepath=CNDO_FILEPATH, box_size=150.)


def update_dashboard(dash_app):
    # First we update the dashboard with the tabs of data we want to show the user:
    dash_app.layout = html.Div([
        dcc.Link(
            html.Button('Return to Homepage', id='home-button', n_clicks=0),
            href='/home',
            refresh=True),
        dcc.Tabs(id="tabs-example-graph", value='tab-1-example-graph', children=[
            dcc.Tab(label='2D Pareto Plots', value='pareto-vs-shape'),
            dcc.Tab(label='Generative Animation Video', value='optimization-search-animation'),
            dcc.Tab(label='Other MOSA Output Figures', value='MOSA-output'),
        ]),
        html.Div(id='tabs-content-example-graph'),
        dcc.Store(id='opt-data-store', storage_type='memory'),
    ])

    @callback(Output('tabs-content-example-graph', 'children'),
              [Input('tabs-example-graph', 'value'), Input('opt-data-store', 'data')])
    def render_content(tab, stored_data):
        opt_data = return_read_data()
        if tab == 'pareto-vs-shape':
            max_mins = find_min_and_maxes(l1=opt_data.archive_post_temperature_initialization,
                                          l2=list(opt_data.MOSA_archive.keys()))

            rep_vs_buck, rep_vs_faces, buck_vs_faces = opt_data.paretos_for_dash_app(
                data_points=list(opt_data.MOSA_archive.keys()),
                max_mins=max_mins)

            return html.Div([
                html.P(
                    "Below are plots of the final Pareto Front found during the MOSA generative process. The red points in each plot specify the 'Pareto Optimal' points for that chart, and the gray points are other data points stored in the archive."),
                html.P(
                    "You are able to left click any of the points in the 2D charts and the 3D visualization for the wireframe structure will appear on the right hand side. A green star represents that currently displayed 3D figure (if a point is left clicked), and the star will appear on all three charts so you know where a design lies across all three paretos."),
                html.P(
                    "The export button (which only works if a point was selected) will export out a PLY representation for the wireframe, as well as oxDNA .top and .dat files."),
                html.P(
                    "You may copy and paste a custom scaffold sequence before exporting a design. Note that if the input design is not long enough for a generated design, then M13 will be used (if shorter than 7,249 bps) or a lambda phage will be used if M13 is too short (by default from DAEDALUS). This should have been controlled during the generation of the input file."),
                html.Div(id='stacked_plots', children=[
                    html.Div(id='rep__vs__buck', children=[dcc.Graph(figure=rep_vs_buck, id='rep_vs_buck')],
                             style={'width': '100%'}),
                    html.Div(id='rep__vs__faces', children=[dcc.Graph(figure=rep_vs_faces, id='rep_vs_faces')],
                             style={'width': '100%'}),
                    html.Div(id='buck__vs__faces', children=[dcc.Graph(figure=buck_vs_faces, id='buck_vs_faces')],
                             style={'width': '100%'})
                ], style={'width': '50%', 'display': 'inline-block'}),

                html.Div(id='pareto_front', style={'width': '35%', 'display': 'inline-block', 'float': 'right'})
            ])

        elif tab == 'optimization-search-animation':
            # all_figs = opt_data.create_pareto_animation(save_or_show='save')
            all_figs = opt_data.all_figures_for_animation

            return html.Div([
                html.H3('Pareto Animation (archive of solutions during optimization process)'),
                html.P(
                    'This animation shows how the design space was searched, how non-dominating solutions are found, and how the trade off (or Pareto Front) was optimized during the MOSA process'),
                html.Div([dcc.Graph(id='animation_fig', figure=all_figs[0], style={'height': '75vh'},
                                    config={'displayModeBar': False})]),
            ])

        elif tab == 'MOSA-output':
            fig1 = opt_data.plot_acceptance_probabilities(save_or_show='save')
            fig2 = opt_data.plot_temperature_profiles(save_or_show='save')
            return html.Div([
                html.H3('Worse Move Acceptance Probability during Generative Process:'),
                html.Div([dcc.Graph(id='animation_fig', figure=fig1)]),

                html.H3('Individual Objective Function Temperatures during Generative Process:'),
                html.Div([dcc.Graph(id='animation_fig', figure=fig2)]),
            ])

    @callback(
        [Output(component_id='pareto_front', component_property='children'),
         Output(component_id='rep__vs__buck', component_property='children'),
         Output(component_id='rep__vs__faces', component_property='children'),
         Output(component_id='buck__vs__faces', component_property='children')],
        # Update the 'ply-rep' figure with the returned design figure
        [Input(component_id='rep_vs_buck', component_property='clickData'),
         Input(component_id='rep_vs_faces', component_property='clickData'),
         Input(component_id='buck_vs_faces', component_property='clickData'),
         Input('opt-data-store', 'data')]
    )
    def update_output_plot(clickData_buck_rep, clickData_rep_faces, clickData_buck_faces, data_file):
        # Check which input triggered the callback
        opt_data = return_read_data()
        triggered_id = callback_context.triggered[0]['prop_id'].split('.')[0]

        if triggered_id == 'rep_vs_buck':
            clickData = clickData_buck_rep
        elif triggered_id == 'rep_vs_faces':
            clickData = clickData_rep_faces
        elif triggered_id == 'buck_vs_faces':
            clickData = clickData_buck_faces
        else:
            clickData = None

        if clickData is not None:
            selected_point_info = clickData['points'][0]
            # First extract the x and y value:
            x, y = selected_point_info['x'], selected_point_info['y']
            # Next, we need to find (from the saved archive), which tuple this point correlates to:
            save_key = ''
            if triggered_id == 'rep_vs_buck':
                for key, value in opt_data.MOSA_archive.items():
                    if x == key[0] and y == key[1]:
                        found_design = value  # Once we find this we save the design and break the loop
                        save_key = key
                        break
            elif triggered_id == 'rep_vs_faces':
                for key, value in opt_data.MOSA_archive.items():
                    if x == key[2] and y == key[1]:
                        found_design = value  # Once we find this we save the design and break the loop
                        save_key = key
                        break
            elif triggered_id == 'buck_vs_faces':
                for key, value in opt_data.MOSA_archive.items():
                    if x == key[2] and y == key[0]:
                        found_design = value  # Once we find this we save the design and break the loop
                        save_key = key
                        break
            # Store the found design in a global variable and calculate the related objective function values:
            session['active_design_to_store'] = save_key
            buck, rep, fac = opt_data.objective_function(design_space=found_design)

            # Next, with the found design, we simply obtain the 3D figure plot to update the ply-rep:
            design_figure = found_design.display_graph(return_figure=True, use_in_dash=True)

            design_figure2 = found_design.display_graph_as_cylindrical_rep(return_figure=True)

            # Now update the selected point plot:
            max_mins = find_min_and_maxes(l1=opt_data.archive_post_temperature_initialization,
                                          l2=list(opt_data.MOSA_archive.keys()))

            design_figure.update_layout(
                margin=dict(l=0, r=75, b=0, t=0),
                width=700,
                height=650
            )
            design_figure2.update_layout(
                margin=dict(l=0, r=75, b=20, t=0),
                width=700,
                height=650
            )

            # Retrieve the chart data again and add in the selected point to the plots:
            rep_vs_buck, rep_vs_faces, buck_vs_faces = opt_data.paretos_for_dash_app(
                data_points=list(opt_data.MOSA_archive.keys()),
                max_mins=max_mins)
            rep_vs_buck.add_trace(go.Scatter(x=[buck], y=[rep], mode='markers',
                                             marker=dict(color='green', line_width=1, size=14, symbol='star'),
                                             opacity=1.0, showlegend=False))
            rep_vs_faces.add_trace(go.Scatter(x=[fac], y=[rep], mode='markers',
                                              marker=dict(color='green', line_width=1, size=14, symbol='star'),
                                              opacity=1.0, showlegend=False))
            buck_vs_faces.add_trace(go.Scatter(x=[fac], y=[buck], mode='markers',
                                               marker=dict(color='green', line_width=1, size=14, symbol='star'),
                                               opacity=1.0, showlegend=False))
            return [dcc.Graph(id='ply_rep', figure=design_figure, config={'displayModeBar': False}),
                    dcc.Graph(id='cyl_rep', figure=design_figure2, config={'displayModeBar': False}),
                    dcc.Input(id='dna-sequence-input', type='text', placeholder='Custom Scaffold Sequence (optional)',
                              style={'margin': 'auto', 'width': '100%'}),
                    dcc.Store(id='dna-sequence-store'),
                    html.Button('Export Design', id='export-design', n_clicks=0,
                                style={'margin': 'auto', 'display': 'block'}),
                    dcc.Download(id="download")], [
                       dcc.Graph(id='rep_vs_buck', figure=rep_vs_buck)], [
                       dcc.Graph(figure=rep_vs_faces, id='rep_vs_faces')], [
                       dcc.Graph(figure=buck_vs_faces, id='buck_vs_faces')]

        else:
            # If nothing is clicked we just show the original plots:
            max_mins = find_min_and_maxes(l1=opt_data.archive_post_temperature_initialization,
                                          l2=list(opt_data.MOSA_archive.keys()))
            rep_vs_buck, rep_vs_faces, buck_vs_faces = opt_data.paretos_for_dash_app(
                data_points=list(opt_data.MOSA_archive.keys()),
                max_mins=max_mins)
            return [], [dcc.Graph(figure=rep_vs_buck, id='rep_vs_buck')], [
                dcc.Graph(figure=rep_vs_faces, id='rep_vs_faces')], \
                   [dcc.Graph(figure=buck_vs_faces, id='buck_vs_faces')]

    @callback(
        Output("download", "data"),
        [Input("export-design", "n_clicks"), Input("dna-sequence-store", "data")]
    )
    def generate_and_download_file(n_clicks, inputSequence):
        if n_clicks is not None and n_clicks > 0:
            # Create a temporary directory:
            savepath = os.path.join(os.getcwd(), 'temp_stored_data')  # Top level save path folder
            random_string = ''.join(random.choice(string.ascii_letters + string.digits) for _ in range(18))
            random_string2 = 'MOSA_Output_Files'
            temp_savepath = os.path.join(savepath, random_string)
            actual_savepath = os.path.join(temp_savepath, random_string2)
            os.makedirs(actual_savepath, exist_ok=True)  # We will always have to create this

            # Dictionary for re-writing folderpath names:
            rename_dict = {random_string: 'MOSA_Output_Files_Zipped'}
            opt_data = return_read_data()
            active_design_key = session['active_design_to_store']
            active_design_to_store = opt_data.MOSA_archive[active_design_key]

            # Save out PLY file:
            writepath = os.path.join(actual_savepath, 'GeneratedDesign_CAD_File.ply')
            active_design_to_store.write_ply_file(filepath=writepath)
            # Save out caDNAno and oxDNA files:
            shortest_edge = active_design_to_store.print_shortest_edge_length()
            # Function call of cando / cadnano export:
            write_out_cando_cadnano(filepath=actual_savepath, plypath=writepath, edge_length=shortest_edge,
                                    dna_seq=inputSequence)

            # Replace string names and zip the file for export:
            shutil.make_archive(temp_savepath.replace(random_string, rename_dict[random_string]), "zip", temp_savepath)

            # Read the ZIP file as bytes
            filezipname = rename_dict[random_string] + '.zip'
            fileZipPath = os.path.join(savepath, filezipname)

            # Create send back and delete temporary folders:
            return_file = dcc.send_file(fileZipPath)
            shutil.rmtree(temp_savepath)
            os.remove(fileZipPath)
            return return_file

        # If the button hasn't been clicked we just return None
        return None

    @callback(
        Output('dna-sequence-store', 'data'),
        Input('dna-sequence-input', 'value'),
        State('dna-sequence-store', 'data')
    )
    def save_text(text_input, stored_data):
        if text_input:
            stored_data = text_input
        else:
            stored_data = "NONE"
        # JSONify:
        return json.dumps(stored_data)

    @callback(
        Output('home-button', 'n_clicks'),
        Input('home-button', 'n_clicks')
    )
    def delete_file(n_clicks):
        # This callback will automatically delete any file stored on the server to make sure we keep "open space"
        if n_clicks is not None and n_clicks > 0:
            filename = session['data_filepath']
            filepath = os.path.join(current_app.config['UPLOAD_FOLDER'], filename)
            if os.path.exists(filepath):
                os.remove(filepath)
            return 0
        else:
            return n_clicks






