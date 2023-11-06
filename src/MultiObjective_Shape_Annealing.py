"""
AJ Vetturini
IDIG and MMBL
Carnegie Mellon University
Advised By: Jon Cagan and Rebecca Taylor

This script incorporates the MOSA algorithm as developed by Suppapitnam, Seffen, Parks, and Clarkson in "A Simulated
Annealing Algorithm for Multiobjective Optimization" (Engineering Optimization, 2000, Vol 33. pp. 59-85)

There are some modifications made to the algorithm in reference to Suppapitnarm, Parks, Shea, and Clarkson in
"Conceptual Design of Bicycle Frames by Multiobjective Shape Annealing"
"""
# Import Modules:
import copy
from random import seed, random, choice, uniform
import warnings
import pickle
import src.Ply_Graphs_3D as pg3
import numpy as np
from dataclasses import dataclass, field
import os
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import pandas as pd
import time
import ipywidgets as widgets
from IPython.display import display

@dataclass
class ShapeAnnealer(object):
    NP_regions: list
    max_x: float
    max_y: float
    max_z: float
    binding_regions: list

    # Parameters that have default values:
    cooling_schedule: str = 'Huang'  # Default to the geometric cooling schedule
    acceptance_function: str = 'standard'  # Values are 'logistic' 'standard' and 'linear'
    routing_software: str = 'DAEDALUS'
    min_edge_length: int = 42
    extend_rule_distance: float = 1.5
    max_number_bps: int = 7249
    max_edge_multiplier: float = 8
    min_face_angle: float = 8.
    check_scaffold_length_constraint: bool = True
    allow_self_intersecting_geometry: bool = False
    max_number_edges_at_node: int = 6
    font_size_dict = {'save': 32, 'show': 32}
    repulsion_distance: float = 2.
    numDecimals: int = 5
    max_number_of_epochs: int = 50

    ## Dummy variables:
    T_buckling: float = 100.  # These are NOT controlled by the user
    T_repulsion: float = 100.  # They are dummy values that are referenced across functions
    T_scaffold: float = 100.  # They are initialized here so I can access them across the data class more easily.

    ## Hyperparameters that can be tuned by user:
    NT1: int = 2000  # (1000 recommended for Initial temperature calculation as recommended by authors of MOSA)
    NT2: int = 500
    Na: int = 200  # Inner loop control, = 0.4NT2 per MOSA algo
    r_b: float = 0.9  # Return to base parameter in (0, 1) and dictates frequency of return (close to 1 == more returns)
    N_Bi: int = 1000  # 2 * NT2 per Shea paper at top. Controls initial # of returns. Larger value == deeper search
    N_Bi_LowerLimit: int = 10  # We need a lower limit for N_Bi which is recommended to be 10 by Shea
    minimal_candidate_set_size: int = 10
    r_i: float = 0.95  # Return fraction parameter, recommended at least 0.9 but can be in (0, 1)
    max_time_of_simulation_minutes: int = 60  # Number of minutes before exiting the optimization for cost protection
    #use_number_faces_as_random_float: bool = True
    C_repulsion: float = 2.0  # Passed in to the DesignSpace, this is default value


    # Cooling Schedule Parameters depending on schedule used (Huang is recommended)
    cooling_rate_geometric: float = 0.8  # The cooling schedule for geometric (only used if Geoemtric schedule is used)
    delta_T: float = 0.8  # Used with the Triki annealing schedule and controls adaptation. May need to modify this.

    T_min: float = 0.00000001  # Because of math the temperature will never truly be 0 so we need a stop point
    random_seed: int = 500  # RNG Seed

    # Archive dictionary:
    MOSA_archive: dict = field(default_factory=dict)
    archive_post_temperature_initialization: list = field(default_factory=list)
    pareto_animation: dict = field(default_factory=dict)

    ## Variables the user should NEVER change that are used for analysis of Shape Annealing:
    rules_selected_tracker: dict = field(default_factory=dict)
    rules_passed_tracker: dict = field(default_factory=dict)
    rules_accepted_tracker: dict = field(default_factory=dict)
    phi_r: float = 1  # Parameter that will be used across algorithm for return to base calcs, do not change from 1!

    # Plotting lists:
    objectives: list = field(default_factory=list)  # List to store all objective values for plotting
    buckling_obj_at_Ti: list = field(default_factory=list)
    repulsion_obj_at_Ti: list = field(default_factory=list)
    Nf_obj_at_Ti: list = field(default_factory=list)
    obj_buckling: list = field(default_factory=list)
    obj_repulsion: list = field(default_factory=list)
    obj_scaffold: list = field(default_factory=list)
    temp_buckling: list = field(default_factory=list)
    temp_repulsion: list = field(default_factory=list)
    temp_scaffold: list = field(default_factory=list)
    acceptance_tracker: list = field(default_factory=list)
    all_figures_for_animation: list = field(default_factory=list)
    total_checked: int = 0
    accepted_via_probability: int = 0

    p_accept_list: list = field(default_factory=list)  # Used to track acceptance probability thru Metropolis

    final_time = None  # Used to just mark the final time in the data object for plotting later

    @staticmethod
    def turn_on_notifier():
        pass

    @staticmethod
    def initialize_shape(design_space: pg3.DesignSpace3D) -> None:
        """
        This will be in charge of calling the initialized shape for Shape Annealing.
        :return: Nothing, it just stores to our data object
        """
        design_space.generate_start_shape()
        design_space.verify_initial_geometry()


    @staticmethod
    def pick_rule() -> str:
        # Rule List:
        rule_list = ['Divide Face', 'Merge Face', 'Retriangulate Face', 'Extend Vertex', 'Normal Direction']
        return choice(rule_list)

    @staticmethod
    def archive_datapoint(archive: list, test_point: tuple) -> tuple[bool, list]:
        """
        This function will take in the current archived points as well as a test point to determine if I need to archive
        the point (and do this archiving if needed).
        """
        # Internal function to determine if a point is dominating or not:
        def dominates(new_point, archived_point):
            # p1 dominates p2 if it's better in at least one objective and not worse in any objective
            # unpacking but in the future I should make this code more efficient, right now i have it for readability
            f1_new, f2_new, f3_new = new_point
            f1_old, f2_old, f3_old = archived_point
            # First line makes sure the new_point is at least as good as archived point
            # Second line check if it IS BETTER than any of f2
            if (f1_new <= f1_old and f2_new <= f2_old and f3_new <= f3_old) and \
                    (f1_new < f1_old or f2_new < f2_old or f3_new < f3_old):
                return True  # If true then this point dominates another!
            else:
                return False

        ### CASE 1: If a candidate solution dominates any member of the archive, those members are removed and the new
        #           solution is added.
        # Find dominated points in the archive and mark their indices for removal:
        to_remove = []
        for i, archive_point in enumerate(archive):
            if dominates(test_point, archive_point):
                # This case means that the test point dominates the archive point meaning we need to remove archive
                # point from the pareto:
                to_remove.append(archive_point)
        # If to_remove has any points added then we know that our archive point dominates at least 1 point and is
        # therefore added to the archive (by returning True) and we remove the points in the to_remove list.
        if to_remove:
            return True, to_remove

        ### CASE 2: If a new solution is dominated by any members of the archive, it is not archived:
        for i, archive_point in enumerate(archive):
            if dominates(archive_point, test_point):
                # If any solution in the archive is dominating the test_point, we return False to not add to archive.
                return False, to_remove

        ### CASE 3: If a test_point does not dominate any solution in the archive, and it is not dominated by any
        #           member of the archive, we DO store the point, but we do NOT remove any from the list. However, here
        #           the to_remove will ALWAYS by of value [] (empty) so we can just return it like that:
        return True, to_remove


    @staticmethod
    def calculate_pareto_front(df: pd.DataFrame, field1: str, field2: str) -> pd.DataFrame:
        """
        This function will take in the pareto points we are plotting and return a dataframe of the optimal points so
        that the pareto front can be plotted. This only calculates the pareto for a 2D set of functions!
        """
        # Initialize a boolean mask to identify non-dominated points
        is_non_dominated = np.ones(len(df), dtype=bool)

        for i, row in df.iterrows():
            if is_non_dominated[i]:
                # Check if there is any other point that dominates the current point
                for j, other_row in df.iterrows():
                    if i != j:  # Skip comparing the point with itself
                        # For two points, A = (x1, y1) B = (x2, y2) iff x1 <= x2 and y1<=y2 then B dominates A
                        #if row[field1] <= other_row[field1] and row[field2] <= other_row[field2]:
                        if (row[field1] <= other_row[field1] and row[field2] <= other_row[field2] and (
                                row[field1] < other_row[field1] or row[field2] < other_row[field2])):
                            is_non_dominated[i] = False  # Mark it as dominated
                            break  # If a point is dominated, break from this inner loop

        # Use the mask to filter the dataframe for non-dominated points
        non_dominated_df = df[is_non_dominated]
        return non_dominated_df


    @staticmethod
    def calculate_3D_pareto(df: pd.DataFrame) -> pd.DataFrame:
        """
        This function calculates the 3D pareto front using the concept of dominance
        """
        pareto_front = []
        # Convert the DataFrame to a NumPy array for faster calculations
        data = df[['Repulsion', 'Buckling', 'Scaffold']].to_numpy()
        # Get the number of data points and dimensions
        num_points, num_dimensions = data.shape
        for i in range(num_points):
            is_pareto = True
            for j in range(num_points):
                if i != j:
                    # Check if point i dominates point j
                    if np.all(data[i] <= data[j]) and np.any(data[i] < data[j]):
                        is_pareto = False
                        break

            if is_pareto:
                pareto_front.append(data[i])

        # Create a new DataFrame with the Pareto front points
        pareto_df = pd.DataFrame(pareto_front, columns=['Repulsion', 'Buckling', 'Scaffold'])

        return pareto_df


    def objective_function(self, design_space: pg3.DesignSpace3D) -> tuple[float, float, float]:
        """
        This function calculates the objective function (or cost function) value and returns it for plotting purposes
        as well as acceptance criterion.
        :return: Float value of the objective function
        """
        total_buckling, F_rep, N_basepairs, N_e = design_space.minimize_RBD_forces()
        avg_buckling = total_buckling / N_e  # Divide the total buckling by number of edges for an average value
        # Since we do not have a composite objective function we do not worry about weights of function values so
        # we return the avg buckling force, the total repulsive force, and 1 / Scaffold usage
        # We use a constant of 100 / Nf strictly to scale temperature values so they are not small in terms of magnitude

        # This term should drive the structure to use "as much" of the (estimated) scaffold as possible up to the limit
        nBs = 1 / (N_basepairs / self.max_number_bps)  # We want to min 1 / (N_i / N_max))

        return round(avg_buckling,self.numDecimals), round(F_rep,self.numDecimals), round(nBs, self.numDecimals)


    def p_accept(self, new: float, old: float, temperature: float) -> float:
        """
        This function will calculate the individual probability for MOSA
        """
        # Use this to catch computational errors in the exponential:
        if self.acceptance_function == 'standard' or self.acceptance_function == 'Standard':
            # I first check to see if new - old / temperature is greater than a threshold as we are passing these to a
            # exponential function. This is because I want to avoid the numpy "very large value" warnings as well as
            # because i want to plot these probability values. A very large value for temp means that the move we are
            # evaluating is VERY positive to the objective function we are testing (it's temperature).

            if ((new - old) / temperature) <= 4.6:
                # This threshold of 4.6 correspond to about ~100 (e^-(new-old/temperature) = e^-4. = 100)
                # which i will return instead for plotting and computatoinal reasons
                return 100
            else:
                return np.exp(((-1 * (new - old)) / temperature))

        elif self.acceptance_function == 'logistic' or self.acceptance_function == 'Logistic':
            if ((new - old) / temperature) <= 4.6:
                # This threshold of 4.6 correspond to very large values of e, and the larger this e value is the actual
                # value of temp used in Logistic APF goes to 0 since it uses e(x) instead of e(-x)
                temp = 0
            else:
                temp = np.exp((new - old) / temperature)  # Logistic Curve Acceptance Probabilty Function
            return 2 / (1 + temp)

        elif self.acceptance_function == 'linear' or self.acceptance_function == 'Linear':
            if ((new - old) / temperature) <= 4.6:
                return 1  # If the threshold is met, we return the value 1 since the APF is min([1, exp])
            else:
                temp = np.exp(((-1 * (new - old)) / temperature))  # Scalar linear Acceptance Probabilty Function
                return min([1, temp])
        else:
            raise Exception('Invalid acceptance probability function input')


    def probability_acceptance(self, old_design_iter: pg3.DesignSpace3D, new_design_iter: pg3.DesignSpace3D) -> float:
        """
        This function compares to design states and calculates the product probability of acceptance to archive a design
        stage per the MOSA implementation
        """
        # Find the objective function values of both design stages:
        old_buckling, old_repulsion, old_scaffold = self.objective_function(old_design_iter)
        new_buckling, new_repulsion, new_scaffold = self.objective_function(new_design_iter)

        # Calculate probabilities:
        p1 = self.p_accept(new=new_buckling, old=old_buckling, temperature=self.T_buckling)
        p2 = self.p_accept(new=new_repulsion, old=old_repulsion, temperature=self.T_repulsion)
        p3 = self.p_accept(new=new_scaffold, old=old_scaffold, temperature=self.T_scaffold)
        final_probability = round(p1 * p2 * p3, self.numDecimals)
        return final_probability


    def quench_temperature(self) -> None:
        """
        This function is in charge of cooling the temperature in MOSA with a few different cooling schedules as a
        potential use
        :return:
        """
        if self.cooling_schedule == 'geometric' or self.cooling_schedule == 'Geoemtric':
            self.T_scaffold *= self.cooling_rate_geometric
            self.T_buckling *= self.cooling_rate_geometric
            self.T_repulsion *= self.cooling_rate_geometric

        elif self.cooling_schedule == 'triki' or self.cooling_schedule == 'Triki':
            # Using a worst-case-scenario where all measured objectives are the same value, use a constant alpha of .8
            if all(x == self.Nf_obj_at_Ti[0] for x in self.Nf_obj_at_Ti):
                alpha_scaffold = 0.5
            else:
                alpha_scaffold = (1 - ((self.T_scaffold * self.delta_T) / np.std(self.Nf_obj_at_Ti)**2))
                # Safe checking:
                if alpha_scaffold < 0.5:
                    alpha_scaffold = 0.5

            if all(x == self.buckling_obj_at_Ti[0] for x in self.buckling_obj_at_Ti):
                alpha_buckling = 0.5
            else:
                alpha_buckling = (1 - ((self.T_buckling * self.delta_T) / np.std(self.buckling_obj_at_Ti) ** 2))
                # Safe checking:
                if alpha_buckling < 0.5:
                    alpha_buckling = 0.5


            if all(x == self.repulsion_obj_at_Ti[0] for x in self.repulsion_obj_at_Ti):
                alpha_repulsion = 0.5
            else:
                alpha_repulsion = (1 - ((self.T_repulsion * self.delta_T) / np.std(self.repulsion_obj_at_Ti) ** 2))
                # Safe checking:
                if alpha_repulsion < 0.5:
                    alpha_repulsion = 0.5

            self.T_scaffold = self.T_scaffold * alpha_scaffold
            self.T_buckling = self.T_buckling * alpha_buckling
            self.T_repulsion = self.T_repulsion * alpha_repulsion
            # We also need to reset the lists after cooling:
            self.Nf_obj_at_Ti, self.buckling_obj_at_Ti, self.repulsion_obj_at_Ti = [], [], []

        elif self.cooling_schedule == 'Huang' or self.cooling_schedule == 'huang':
            if all(x == self.Nf_obj_at_Ti[0] for x in self.Nf_obj_at_Ti):
                alpha_scaffold = 0.5
            else:
                alpha_scaffold = max((0.5, np.exp(-0.7 * self.T_scaffold / np.std(self.Nf_obj_at_Ti))))

            if all(x == self.buckling_obj_at_Ti[0] for x in self.buckling_obj_at_Ti):
                alpha_buckling = 0.5
            else:
                alpha_buckling = max((0.5, np.exp(-0.7 * self.T_buckling / np.std(self.buckling_obj_at_Ti))))

            if all(x == self.repulsion_obj_at_Ti[0] for x in self.repulsion_obj_at_Ti):
                alpha_repulsion = 0.5
            else:
                alpha_repulsion = max((0.5, np.exp(-0.7 * self.T_repulsion / np.std(self.repulsion_obj_at_Ti))))

            # After calculating the proper alpha we decrease temperature:
            self.T_scaffold *= alpha_scaffold
            self.T_buckling *= alpha_buckling
            self.T_repulsion *= alpha_repulsion

            # Finally, since temperature has now cooled, we reset the lists for the standard deviation calcs:
            self.Nf_obj_at_Ti, self.buckling_obj_at_Ti, self.repulsion_obj_at_Ti = [], [], []

        else:
            raise Exception('Invalid cooling schedule selected, options are: geometric, triki, huang')

        # Update temperature lists for plotting:
        self.temp_buckling.append(self.T_buckling)
        self.temp_repulsion.append(self.T_repulsion)
        self.temp_scaffold.append(self.T_scaffold)


    def calculate_and_store_objective_functions(self, design_space: pg3.DesignSpace3D) -> None:
        """
        This function just stores the objective function values since I do this often enough and want to keep the code
        efficient
        """
        buckling, repulsion, scaffold = self.objective_function(design_space=design_space)
        # Store these values AND archive the initial state!
        self.buckling_obj_at_Ti.append(buckling)
        self.repulsion_obj_at_Ti.append(repulsion)
        self.Nf_obj_at_Ti.append(scaffold)
        archive_data, remove_list = self.archive_datapoint(archive=list(self.MOSA_archive.keys()),
                                                           test_point=(buckling, repulsion, scaffold))
        if archive_data:
            # Archive is stored as (magnitudes of objective functions) : design_space object
            self.MOSA_archive[(buckling, repulsion, scaffold)] = copy.deepcopy(design_space)

            # This removal loop will only be needed if we are archiving a point, if we are not archiving a point we will
            # never be removing a point
            for remove_from_archive in remove_list:
                self.MOSA_archive.pop(remove_from_archive)


    def calc_and_store_obj_func_values(self, design_space: pg3.DesignSpace3D) -> None:
        buckling, repulsion, scaffold = self.objective_function(design_space=design_space)
        # Store these values AND archive the initial state!
        self.obj_buckling.append(buckling)
        self.obj_repulsion.append(repulsion)
        self.obj_scaffold.append(scaffold)


    def initialize_temperatures(self, design_space: pg3.DesignSpace3D):
        """
        This function is used to start the MOSA algorithm and will initially search the indvidiual design spaces to find
        trade offs.
        """
        # Start by evaluating the current start shape:
        self.calculate_and_store_objective_functions(design_space=design_space)

        # We will ALWAYS accept a move, but we need a counter and while loop since we can still violate constraints
        accepted_moves = 0
        while_loop_counter = 0
        Nb_reset = 100  # Test value for return to base on the initial walk through space for a better start space.
        while accepted_moves < self.NT1:
            # We first randomly apply a move and verify it doesn't violate any design constraints:
            production_rule_choice = self.pick_rule()
            next_iter_design_space = copy.deepcopy(design_space)
            next_iter_design_space.production_rules(rule_selection=production_rule_choice)
            if next_iter_design_space.master_design_constraint_checker():
                # If we pass the design constraints we ALWAYS accept as we are just walking the design space randomly
                # to get a starting temperature based on literature
                self.calculate_and_store_objective_functions(design_space=next_iter_design_space)
                accepted_moves += 1  # Add one to accepted moves as we map out the design space.
                if accepted_moves % Nb_reset == 0:
                    # When we hit threshold for the reset on the initial walk through space, we reset to a random state
                    # to continue generation
                    X_N1 = self.return_to_base(initial_search=True)
                    del design_space
                    design_space = copy.deepcopy(X_N1)
                else:
                    del design_space
                    design_space = copy.deepcopy(next_iter_design_space)  # We also reassign the design space

            # For the sake of no infinite loops, I use a while loop counter here to ensure it doesn't get stuck:
            while_loop_counter += 1
            if while_loop_counter > 100000:
                warnings.warn(f'Unable to apply {self.NT1} without violating design constraints which may lead to '
                              f'an inferior solution, be aware!')
                accepted_moves = self.NT1 + 1  # Just break the loop condition...

        # After the while loop has terminated, we use the methodology by S. R. White from "Concepts of scale in
        # simulated annealing. IEEE International Conference of Computer Aided Design, Port Chester, New York.
        # pp. 646-651" to initialize each temperature of each objective function:
        self.T_buckling = round(float(np.std(self.buckling_obj_at_Ti)), self.numDecimals)
        self.T_repulsion = round(float(np.std(self.repulsion_obj_at_Ti)), self.numDecimals)
        self.T_scaffold = round(float(np.std(self.Nf_obj_at_Ti)), self.numDecimals)

        # We append the objectives to the master list for plotting and reset these lists used in these NT1 steps
        self.obj_buckling.extend(self.buckling_obj_at_Ti)
        self.obj_repulsion.extend(self.repulsion_obj_at_Ti)
        self.obj_scaffold.extend(self.Nf_obj_at_Ti)
        self.Nf_obj_at_Ti, self.buckling_obj_at_Ti, self.repulsion_obj_at_Ti = [], [], []

        # We will plot the pareto front of this process at the end to see how the design evolves:
        self.archive_post_temperature_initialization = list(self.MOSA_archive.keys())


    def keep_annealing(self, curTime, numEpochs) -> bool:
        """
        This function will be used to control the end conditions for the Simulated Annealing loop
        """
        elapsed_time = time.time() - curTime
        hours, remainder = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        time_formatted = "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)

        ### First check if temperatures are low enough to end the run:
        if self.T_scaffold < self.T_min and self.T_repulsion < self.T_min and self.T_buckling < self.T_min:
            # If all temperatures drop below the minimal threshold, we end the search
            print(f'MOSA ended due to temperatures dropping below threshold of {self.T_min}')
            self.final_time = time_formatted
            return False
        ### Next check if max time has elapsed
        elif (time.time() - curTime) >= (self.max_time_of_simulation_minutes * 60):
            # If the current elapsed time is greater than the threshold in seconds (why i multiply by 60 above), we end
            print(f'MOSA ended due to max time input of {self.max_time_of_simulation_minutes} minutes')
            self.final_time = time_formatted
            return False
        elif numEpochs > self.max_number_of_epochs:
            print(f'MOSA ended due to max number of epochs of {self.max_number_of_epochs} being met')
            self.final_time = time_formatted
            return False
        else:
            # Otherwise we return True meaning we keep looping
            return True


    def calculate_pareto_isolation(self, initial_search: bool = False) -> pd.DataFrame:
        """
        This function is used to rank and calculate the isolations of the pareto front during optimization
        """
        pareto = list(self.MOSA_archive.keys())
        # First I am going to split the pareto into lists for each objective function:
        buckling, repulsion, scaffold = [list(col) for col in zip(*pareto)]
        maxBRF = [max(buckling), max(repulsion), max(scaffold)]
        minBRF = [min(buckling), min(repulsion), min(scaffold)]

        # Loop over all solutions
        isolation_dict = {}  # Dictionary to track measures of isolation
        for i in pareto:
            # Our inner loop is for all other solutions in the pareto
            Iij = 0
            for j in pareto:
                if i == j:
                    # Here nothing happens because we don't compare the isolation of the same point
                    continue
                else:
                    # Calculate isolation value for this index:
                    for f1, f2, FMAX, FMIN in zip(i, j, maxBRF, minBRF):
                        if FMAX == FMIN:
                            Iij += 1  # We just add a constant value since repulsion can be fully minimized
                            # The constant will be added to every isolation value since if FMAX and FMIN are the same
                            # then that means all values in i and j for that column are 0.
                        else:
                            Iij += (((f1 - f2) / (FMAX - FMIN))**2)

            # Before continuing, we store this isolation measure:
            isolation_dict[i] = Iij

        # After calculating all isolation values, we then rank in order of decreasing isolation distance WITH THE
        # EXCEPTION of the extreme solutions (i.e. solutions corresponding to extrema in the trade-off / solutions
        # with the lowest objective function. Sorting the dictionary:
        points, isolations = [], []
        for k, v in isolation_dict.items():
            # Since I am currently using a relatively simple objective function space, I will not remove the extrema
            # datapoints, as these do not explicitly constitute a barely feasible solution.
            """if k[0] == minBRF[0] or k[1] == minBRF[1] or k[2] == minBRF[2]:
                # If the point we are checking contains a point in the extrema (which in this case is a MINIMIZATION
                # problem), then we are not going to choose it in our return to base
                continue
            else:
                # If the point is not an extrema, we are going to add it to this 2D list:
                points.append(k)
                isolations.append(v)"""
            points.append(k)
            isolations.append(v)

        # Combining into a dataframe to sort by descending order per MOSA algorithm:
        df = pd.DataFrame({'points': points, 'isolations': isolations})
        df = df.sort_values(by='isolations', ascending=False)
        datapoints = df.shape[0]
        selection_set_size = round(self.phi_r * datapoints)  # Rounding because we need an integer # of points to pick
        # We also incorporate a minimal selection set size to prevent ever focusing in on just a few solutions:
        if selection_set_size < self.minimal_candidate_set_size:
            selection_set_size = self.minimal_candidate_set_size

        # We select the first "selection_set_size" number of points to be into "candidate solutions" we will choose from
        candidate_list = df[:selection_set_size]

        # After creating the candidate list, we update values for future selections. But we only do this when we are not
        # performing the initial walk through space
        if not initial_search:
            self.phi_r *= self.r_i

        return candidate_list


    def return_to_base(self, initial_search: bool = False) -> pg3.DesignSpace3D:
        """
        The return to base functionality is used to expose the trade-off between objective functions. The return to base
        will randomly change the active design state to somewhere else within the Pareto so that we optimize the
        pareto front

        Here we utilize the "intelligent" return-to-base strategy laid out in Conceptual Design of Bicycle Frames by
        Multiobjective Shape Annealing. The purpose is to prefer to "return to base" to extreme solutions (or highly
        isolated solutions) to explore around these unknown areas of the design space
        """
        # First we calculate the paretio isolation parameter to determine where we will "return to base" to search
        candidate_list = self.calculate_pareto_isolation(initial_search=initial_search)

        # Finally, we randomly select a single candidate to become the new active state:
        new_state = candidate_list.sample(n=1)
        new_state = tuple(new_state['points'])  # Convert the datatype...
        new_state_data_object = self.MOSA_archive[(new_state[0][0], new_state[0][1], new_state[0][2])]

        # Lastly, before returning the data object, we update the N_Bi term which dictates how often we will return to
        # base during the search. In general, the longer the search the more often we return to base to exploit the
        # tradeoffs of the objective functions:
        # NOTE: We only lower these values once inside the MOSA loop, for the initial walk we keep the default value.
        if not initial_search:
            # We only lower these values once we start the actual MOSA search
            self.N_Bi = int(round(self.r_b * self.N_Bi))
            if self.N_Bi < self.N_Bi_LowerLimit:
                self.N_Bi = self.N_Bi_LowerLimit  # If we ever go below the lower limit, we just set N_Bi to the limit

        return new_state_data_object


    def attempt_to_archive(self, new_design_space: pg3.DesignSpace3D) -> bool:
        """
        This function will be used to determine if a point is archived or not during the inner While loop of MOSA
        We return True if the point was ACCEPTED and False if the point was REJECTED
        """
        # Calculate objective function values and pass to archive function:
        buckling, repulsion, scaffold = self.objective_function(design_space=new_design_space)
        archive_data, remove_list = self.archive_datapoint(archive=list(self.MOSA_archive.keys()),
                                                           test_point=(buckling, repulsion, scaffold))

        # If we are archiving the data point we use the same logic to update the dictioary depending on dominance
        if archive_data:
            # First store data
            self.buckling_obj_at_Ti.append(buckling)
            self.repulsion_obj_at_Ti.append(repulsion)
            self.Nf_obj_at_Ti.append(scaffold)
            # Archive is stored as (magnitudes of objective functions) : design_space object
            self.MOSA_archive[(buckling, repulsion, scaffold)] = copy.deepcopy(new_design_space)
            # This removal loop will only be needed if we are archiving a point, if we are not archiving a point we will
            # never be removing a point
            for remove_from_archive in remove_list:
                self.MOSA_archive.pop(remove_from_archive)

            # Finally, we return True as if a datapoint is archived, then that state is always accepted per MOSA:
            return True
        # If we do not archive, we return False to do a probability test:
        return False


    def metropolis(self, p_accept: float) -> bool:
        """
        Metropolis-Hastings Algorithm for accepting moves in Shape Annealing.
        :param p_accept: Calculated acceptance probability
        :return: True: Accepted the new move, False: Rejected
        """
        ### PER THE MOSA ALGORITHM: Sometimes, p_accept will be greater than unity (1) for a given function. In this
        #                           case, we will always accept it:
        if p_accept > 1.:
            return True

        # Otherwise, in the case of MINIMIZATION (which all of my objective functions are):
        self.p_accept_list.append(p_accept)  # Append to list for plotting later to verify algorithm
        if random() < p_accept:
            # If the random generated point is less than our acceptance probability, we return True to accept the state
            # Over time, the p_accept should go lower and lower meaning we accept less and less worse states.
            return True
        else:
            # Otherwise we do not accept the state and return False
            return False


    def store_pareto_for_animation(self, epoch: int) -> None:
        """
        This function takes in the current epoch and stores the required data for plotting an animation or video
        to see the growth / exploration of the design space from a pareto POV
        """
        curData = list(self.MOSA_archive.keys())
        self.pareto_animation[epoch] = curData


    def begin_MOSA(self) -> dict:
        """
        This method is responsible for starting MOSA and will contain the code for the MOSA algorithm
        """
        # Turn on the notifier class to make sure AJ is not missing things during development stupidly:
        self.turn_on_notifier()
        # Initialize the design space:
        total_design_space = pg3.DesignSpace3D(maxX=self.max_x, maxY=self.max_y, maxZ=self.max_z,
                                               NP_Regions=self.NP_regions, routing_software=self.routing_software,
                                               binding_regions=self.binding_regions, random_seed=self.random_seed,
                                               min_edge_length=self.min_edge_length,
                                               extend_rule_distance=self.extend_rule_distance,
                                               max_number_basepairs=self.max_number_bps,
                                               max_edge_length_multiplier=self.max_edge_multiplier,
                                               min_face_angle=self.min_face_angle,
                                               check_scaffold_length_constraint=self.check_scaffold_length_constraint,
                                               allow_self_intersecting_geometry=self.allow_self_intersecting_geometry,
                                               max_number_edges_at_node=self.max_number_edges_at_node,
                                               repulsion_distance=self.repulsion_distance)
        # Next we initialize the Random Seeds to be the value used for this Shape Annealing:
        seed(self.random_seed)
        np.random.seed(self.random_seed)
        total_design_space.initialize_random_seed(rng_seed=self.random_seed)
        # Initialize any nonterminal nodes (if any)
        total_design_space.find_nonterminal_nodes()
        # Finally, initialize the geometry we will use for the simulation
        self.initialize_shape(design_space=total_design_space)

        ### STEP 1 of MOSA: Initialize temperatures for all objective functions
        self.initialize_temperatures(design_space=total_design_space)

        ### STEP 2 of MOSA: Use intelligent Return to Base to reset the reset the current design space along the pareto
        X_N1 = self.return_to_base()

        ### STEP 3 of MOSA: Begin loops of Simulated Annealing with temperature full considered above
        inner_count = 0  # Used to track the internal loop count number so we know when to quench
        accepted_count = 0  # Used to track the number of accepted moves for quenching
        return_to_base_count = 0  # Used to track number of loops so we know when to return to base
        curEpoch = 1

        ### BEGIN MOSA (Currently a simple exit condition of Sim. Anneal)
        start_time = time.time()  # We use this to ensure a run doesn't go too long.
        while self.keep_annealing(curTime=start_time, numEpochs=curEpoch):
            # Step 1: Perturb the current design state into a sample state:
            production_rule_choice = self.pick_rule()
            X_N2 = copy.deepcopy(X_N1)
            X_N2.production_rules(rule_selection=production_rule_choice)
            # We always use HARD design constraints in MOSA which changes the flow chart a bit
            if X_N2.master_design_constraint_checker():
                # If we pass the constraints check, we attempt to archive:
                archive_datapoint = self.attempt_to_archive(new_design_space=X_N2)
                if archive_datapoint:
                    del X_N1  # Remove from memory to be safe
                    # IF archived, we accept this X_N2 state as our new state and store the objective function values:
                    X_N1 = copy.deepcopy(X_N2)  # Set X_N1 to X_N2 to signify the acceptance
                    self.calc_and_store_obj_func_values(design_space=X_N1)
                    accepted_count += 1  # If we archive we accept and iterated this counter
                else:
                    # Otherwise we use a probability tet to determine acceptance per MOSA
                    p_accept = self.probability_acceptance(old_design_iter=X_N1, new_design_iter=X_N2)
                    # Now we use metropolis to test the acceptance probability:
                    accept_new_state = self.metropolis(p_accept=p_accept)
                    """
                    Note: Overall, this selection criterion used is a weak point of this MOSA and it is only 
                          recommended for two or three objectives. Further research should be conducted here. 
                    """
                    if accept_new_state:
                        # If we accept via Metropolis then we set the new state to the active state and store objectives
                        del X_N1
                        X_N1 = copy.deepcopy(X_N2)
                        self.calc_and_store_obj_func_values(design_space=X_N1)
                        self.accepted_via_probability += 1
                        self.total_checked += 1
                        accepted_count += 1  # If we probabilistically accept, we iterate accepted_count
                    else:
                        # If we do not accept via Metropolis, we just continue along with X_N1 as the active state
                        self.total_checked += 1
                # Following this, we increment the inner loop count
                inner_count += 1
                if self.total_checked > 1:  # We only start counting once we have passed a constraint check
                    self.acceptance_tracker.append(self.accepted_via_probability / self.total_checked)
            else:
                # If we do not pass the design_constraint check, we still iterate the inner loop count:
                inner_count += 1

            ### PERIODICALLY: Return to base
            return_to_base_count += 1  # Increment the return to base count to determine when we should next return:
            if return_to_base_count == self.N_Bi:
                # If we are returning to base, then we call the new state and set it to XN1:
                returned_state = self.return_to_base()
                X_N1 = returned_state
                return_to_base_count = 0  # Reset this count here

            ### PERIODICALLY: Annealing of Temperature per van Laarhoven and Aarts
            if inner_count >= self.NT2 or accepted_count >= self.Na:
                self.quench_temperature()
                # After quenching we set the counters for BOTH inner_count and accepted_count to 0 for the next temperature
                inner_count, accepted_count = 0, 0
                self.store_pareto_for_animation(epoch=curEpoch)
                curEpoch += 1  # When quenching, we add 1 to epoch.

        # After the while loop concludes / exits, we return the ARCHIVED dataset so we can access any of these values
        return self.MOSA_archive


    def dump_MOSA_object_data(self, savepath: str = None) -> None:
        """
        This function will dump a pickle file to read in later and use for debugging and things
        """
        if savepath is None:
            savepath = os.path.join(os.getcwd(), 'SavedGraph')
            if not os.path.exists(savepath):
                os.mkdir(savepath)

        # I am using the 'aj1' extension because I want to track which versions of my software can generate which
        # types of solutions and plots. For example, if I add in features later on I want to differentiate those types
        # of simulation output.
        output_file = os.path.join(savepath, 'MOSA_export_data.aj1')
        copy_of_data = copy.deepcopy(self)
        with open(output_file, 'wb') as f:
            pickle.dump(copy_of_data, f)


    def plot_pareto_surface_from_archive(self, data_points: list, max_mins: list, save_or_show: str = 'save',
                                         savepath: str = None, plot_title='None') -> None:
        """
        This function will use plotly to plot (and save) the pareto front. We save out the Pareto at the end of the
        initial temperature optimization during White's algorithm and at the end of the MOSA process.
        """
        # We will create side-by-side plots so we can see how the pareto optimal surface came to be during MOSA:
        fig = make_subplots(rows=3, cols=1,
                            subplot_titles=("Repulsion vs Buckling", "Repulsion vs Estimated Scaffold Usage",
                                            "Buckling vs Estimated Scaffold Usage", "3D Pareto Surface"),
                            )


        # Process data to add a singular trace to the figure to make life easier:
        df = pd.DataFrame(data_points, columns=['Buckling', 'Repulsion', 'Scaffold'])
        df_repulsion_buckling = df[['Repulsion', 'Buckling']].copy()
        df_repulsion_scaffold = df[['Repulsion', 'Scaffold']].copy()
        df_buckling_scaffold = df[['Buckling', 'Scaffold']].copy()

        # Find the pareto optimal points:
        pareto_rep_buck = self.calculate_pareto_front(df=df_repulsion_buckling, field1='Buckling', field2='Repulsion')
        pareto_rep_scaffold = self.calculate_pareto_front(df=df_repulsion_scaffold, field1='Repulsion', field2='Scaffold')
        pareto_buck_scaffold = self.calculate_pareto_front(df=df_buckling_scaffold, field1='Buckling', field2='Scaffold')

        #### PLOTTING ALL PARETO POINTS (NOT THE PARETO FRONT):
        fig.add_trace(go.Scatter(x=df['Buckling'], y=df['Repulsion'], mode='markers',
                                 name='Repulsion vs Buckling', marker=dict(color='darkgray', line_width=1, size=10),
                                 opacity=0.6, showlegend=False),
                      row=1, col=1)
        fig.add_trace(go.Scatter(x=df['Scaffold'], y=df['Repulsion'], mode='markers', name='Repulsion vs Estimated Scaffold Usage',
                                 marker=dict(color='darkgray', line_width=1, size=10),
                                 opacity=0.6, showlegend=False),
                      row=2, col=1)
        fig.add_trace(go.Scatter(x=df['Scaffold'], y=df['Buckling'], mode='markers', name='Buckling vs Estimated Scaffold Usage',
                                 marker=dict(color='darkgray', line_width=1, size=10), opacity=0.6,
                                 showlegend=False),
                      row=3, col=1)

        #### PLOTTING ALL PARETO FRONT:
        fig.add_trace(go.Scatter(x=pareto_rep_buck['Buckling'], y=pareto_rep_buck['Repulsion'], mode='markers',
                                 name='Repulsion vs Buckling', marker=dict(color='crimson', line_width=1, size=10),
                                 opacity=1.0, showlegend=False),
                      row=1, col=1)
        fig.add_trace(go.Scatter(x=pareto_rep_scaffold['Scaffold'], y=pareto_rep_scaffold['Repulsion'], mode='markers',
                                 name='Repulsion vs Estimated Scaffold Usage', marker=dict(color='crimson', line_width=1, size=10),
                                 opacity=1.0, showlegend=False),
                      row=2, col=1)
        fig.add_trace(go.Scatter(x=pareto_buck_scaffold['Scaffold'], y=pareto_buck_scaffold['Buckling'], mode='markers',
                                 name='Buckling vs Estimated Scaffold Usage',
                                 marker=dict(color='crimson', line_width=1, size=10), opacity=1.0,
                                 showlegend=False),
                      row=3, col=1)

        fig.update_layout(plot_bgcolor='white')

        count = 1
        for i, j in zip([1, 2, 3], [1, 1, 1]):
            fig.update_xaxes(title_font=dict(size=20), tickfont=dict(size=16), row=i, col=j)
            fig.update_yaxes(title_font=dict(size=20), tickfont=dict(size=16), row=i, col=j)

            subplot_title = fig.layout.annotations[count - 1]['text']
            fig.update_annotations(font=dict(size=24), selector=dict(text=subplot_title))
            count += 1

        # Format of max_mins: [b_max, b_min, r_max, r_min, f_max, f_min]
        # Calculate padding so that points appear on the screen!
        padding_buckling = (max_mins[0] - max_mins[1]) * 0.1
        padding_repulsion = (max_mins[2] - max_mins[3]) * 0.1
        padding_scaffold = (max_mins[4] - max_mins[5]) * 0.1


        fig.update_xaxes(title_text='Buckling', row=1, col=1, linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[1] - padding_buckling, max_mins[0] + padding_buckling])
        fig.update_yaxes(title_text='Repulsion', row=1, col=1, linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[3] - padding_repulsion, max_mins[2] + padding_repulsion])

        fig.update_xaxes(title_text='Scaffold', row=2, col=1, linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[5] - padding_scaffold, max_mins[4] + padding_scaffold])
        fig.update_yaxes(title_text='Repulsion', row=2, col=1, linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[3] - padding_repulsion, max_mins[2] + padding_repulsion])

        fig.update_xaxes(title_text='Scaffold', row=3, col=1, linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[5] - padding_scaffold, max_mins[4] + padding_scaffold])
        fig.update_yaxes(title_text='Buckling', row=3, col=1, linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[1] - padding_buckling, max_mins[0] + padding_buckling])

        if save_or_show == 'show':
            fig.update_layout(title=dict(text='<b>' + plot_title + '</b>', font=dict(size=22, color='blue')))
            fig.show()
        else:
            return fig


    def plot_all_objective_functions(self, save_or_show: str = 'save', savepath: str = None, plot_title='None'):
        """
        This function will use plotly to plot (and save) the pareto front. We save out the Pareto at the end of the
        initial temperature optimization during White's algorithm and at the end of the MOSA process.
        """
        # We will create side-by-side plots so we can see how the pareto optimal surface came to be during MOSA:
        fig = make_subplots(rows=1, cols=3,
                            subplot_titles=("Minimize Avg. Buckling Forces", "Minimize Total Repulsive Forces",
                                            "Minimize 1 / Estimated Scaffold Usage"), start_cell="top-left")
        x = np.arange(1, len(self.obj_buckling) + 1)
        fig.add_trace(go.Scatter(x=x, y=self.obj_buckling, mode='lines',
                                 line=dict(color="black", dash='solid'),
                                 showlegend=False),
                      row=1, col=1)
        fig.add_trace(go.Scatter(x=x, y=self.obj_repulsion, mode='lines',
                                 line=dict(color="black", dash='solid'),
                                 showlegend=False),
                      row=1, col=2)
        fig.add_trace(go.Scatter(x=x, y=self.obj_scaffold, mode='lines',
                                 line=dict(color="black", dash='solid'),
                                 showlegend=False),
                      row=1, col=3)
        fig.update_layout(plot_bgcolor='white')

        count = 1
        for i in [1, 2, 3]:
            fig.update_xaxes(title_font=dict(size=20), tickfont=dict(size=16), row=1, col=i)
            fig.update_yaxes(title_font=dict(size=20), tickfont=dict(size=16), row=1, col=i)

            subplot_title = fig.layout.annotations[count - 1]['text']
            fig.update_annotations(font=dict(size=24), selector=dict(text=subplot_title))
            count += 1

        fig.update_xaxes(title_text='Iteration #', row=1, col=1, linecolor='black', linewidth=1, mirror=True)
        fig.update_yaxes(title_text='Average Buckling Force', row=1, col=1, linecolor='black', linewidth=1, mirror=True)

        fig.update_xaxes(title_text='Iteration #', row=1, col=2, linecolor='black', linewidth=1, mirror=True)
        fig.update_yaxes(title_text='Total Repulsive Forces', row=1, col=2, linecolor='black', linewidth=1, mirror=True)

        fig.update_xaxes(title_text='Iteration #', row=1, col=3, linecolor='black', linewidth=1, mirror=True)
        fig.update_yaxes(title_text='Estimated Scaffold Usage', row=1, col=3, linecolor='black', linewidth=1, mirror=True)
        fig.update_layout(title=dict(text='<b>'+plot_title+'</b>', font=dict(size=22, color='blue')))
        if save_or_show == 'show':
            fig.show()
        else:
            return fig


    def plot_acceptance_probabilities(self, save_or_show: str = 'save', savepath: str=None, plot_title='None') -> None:
        """
        This function will use plotly to plot (and save) the acceptance probabilty during MOSA to verify that we are
        attempting to only accept the objectively good solutions
        """
        # We will create side-by-side plots so we can see how the pareto optimal surface came to be during MOSA:
        fig = make_subplots(rows=1, cols=2,
                            subplot_titles=("Acceptance Probability at Iteration i", "Actual Probability Calculations"),
                            start_cell="top-left")

        x = np.arange(1, len(self.acceptance_tracker) + 1)
        x2 = np.arange(1, len(self.p_accept_list) + 1)
        fig.add_trace(go.Scatter(x=x, y=self.acceptance_tracker, mode='lines',
                                 line=dict(color="black", dash='solid'),
                                 showlegend=False), row=1, col=1)
        fig.add_trace(go.Scatter(x=x2, y=self.p_accept_list, mode='markers',
                                 marker=dict(color='darkgray', line_width=1, size=10), opacity=0.6,
                                 showlegend=False), row=1, col=2)
        fig.add_trace(go.Scatter(x=x2, y=[1] * len(x2),  mode='lines', line=dict(color="red", dash='dash'),
                                 showlegend=False), row=1, col=2)

        fig.update_layout(plot_bgcolor='white')
        # Update axes:
        fig.update_xaxes(title_text='Iteration i', title_font=dict(size=20),
                         tickfont=dict(size=16), linecolor='black', linewidth=1, mirror=True, row=1, col=1)
        fig.update_yaxes(title_text='Acceptance Probability at iteration i', title_font=dict(size=20),
                         tickfont=dict(size=16), linecolor='black', linewidth=1, mirror=True, row=1, col=1)

        fig.update_xaxes(title_text='"Worse Move" Iteration #', title_font=dict(size=20),
                         tickfont=dict(size=16), linecolor='black', linewidth=1, mirror=True, row=1, col=2)
        fig.update_yaxes(title_text='Actual Probability Calculation', title_font=dict(size=20),
                         range=[0, 1.1], tickfont=dict(size=16), linecolor='black', linewidth=1, mirror=True,
                         row=1, col=2)


        if save_or_show == 'show':
            fig.update_layout(title=dict(text='<b>' + plot_title + '</b>', font=dict(size=22, color='blue')))
            fig.show()
        else:
            return fig


    def plot_temperature_profiles(self, save_or_show: str = 'save', savepath: str=None, plot_title='None') -> None:
        """
        This function will use plotly to plot (and save) the objective function temperatures during the run so we can
        verify how the search is going.
        """
        # We will create side-by-side plots so we can see how the pareto optimal surface came to be during MOSA:
        fig = go.Figure()
        x = np.arange(1, len(self.temp_buckling) + 1)
        fig.add_trace(go.Scatter(x=x, y=self.temp_buckling, mode='lines', name='Buckling Temperature',
                                 line=dict(color="black", dash='solid'),
                                 showlegend=True))
        fig.add_trace(go.Scatter(x=x, y=self.temp_repulsion, mode='lines', name='Repulsion Temperature',
                                 line=dict(color="blue", dash='dash'),
                                 showlegend=True))
        fig.add_trace(go.Scatter(x=x, y=self.temp_scaffold, mode='lines', name='Scaffold Temperature',
                                 line=dict(color="green", dash='dot'),
                                 showlegend=True))

        fig.update_layout(plot_bgcolor='white')
        fig.update_xaxes(title_text='Epoch #', title_font=dict(size=20),
                         tickfont=dict(size=16), linecolor='black', linewidth=1, mirror=True)
        fig.update_yaxes(title_text='Temperature at Iteration #', title_font=dict(size=20),
                         tickfont=dict(size=16), linecolor='black', linewidth=1, mirror=True)


        if save_or_show == 'save':
            return fig
        else:
            fig.update_layout(title=dict(text='<b>' + plot_title + '</b>', font=dict(size=22, color='blue')))
            fig.show()


    def create_pareto_animation(self, plot_title='None', save_or_show: str = 'save') -> None:
        """
        This function uses the dictionary to create Frame obejcts to create a plotly animation
        """
        def find_min_and_maxes(l1: list, l2: list) -> list:
            """
            This function finds the min and max values for the repulsion, buckling, and faces between the two lists so that the
            range is FIXED in both plots to understand how MOSA searches the design space
            """
            buckling, repulsion, scaffold = [list(col) for col in zip(*l1)]
            buckling2, repulsion2, scaffold2 = [list(col) for col in zip(*l2)]

            b_max, b_min = max([max(buckling), max(buckling2)]), min([min(buckling), min(buckling2)])
            r_max, r_min = max([max(repulsion), max(repulsion2)]), min([min(repulsion), min(repulsion2)])
            f_max, f_min = max([max(scaffold), max(scaffold2)]), min([min(scaffold), min(scaffold2)])

            return [b_max, b_min, r_max, r_min, f_max, f_min]
        ## CREATE INITIAL FRAME / FIGURE with PRE-MOSA DATA:
        fig = make_subplots(rows=2, cols=2,
                            subplot_titles=("Repulsion vs Buckling", "Repulsion vs Estimated Scaffold Usage",
                                            "Buckling vs Estimated Scaffold Usage", ""),
                            start_cell="top-left",
                            specs=[
                                [{"type": "scatter"}, {"type": "scatter"}],
                                [{"type": "scatter"}, {"type": "scatter"}]
                            ]
                            )
        #fig_3D = go.Figure()  # 3D Pareto
        all_figs = []

        df = pd.DataFrame(self.archive_post_temperature_initialization, columns=['Buckling', 'Repulsion', 'Scaffold'])
        df_repulsion_buckling = df[['Buckling', 'Repulsion']].copy()
        df_repulsion_scaffold = df[['Repulsion', 'Scaffold']].copy()
        df_buckling_scaffold = df[['Buckling', 'Scaffold']].copy()

        # Find the pareto optimal points:
        pareto_rep_buck = self.calculate_pareto_front(df=df_repulsion_buckling, field1='Repulsion', field2='Buckling')
        pareto_rep_scaffold = self.calculate_pareto_front(df=df_repulsion_scaffold, field1='Repulsion', field2='Scaffold')
        pareto_buck_scaffold = self.calculate_pareto_front(df=df_buckling_scaffold, field1='Buckling', field2='Scaffold')

        #### PLOTTING ALL PARETO POINTS (NOT THE PARETO FRONT):
        fig.add_trace(go.Scatter(x=df['Repulsion'], y=df['Buckling'], mode='markers',
                                 name='Buckling vs Repulsion', marker=dict(color='darkgray', line_width=1, size=10),
                                 opacity=0.6, showlegend=False),
                      row=1, col=1)
        fig.add_trace(go.Scatter(x=df['Scaffold'], y=df['Repulsion'], mode='markers', name='Repulsion vs Estimated Scaffold Usage',
                                 marker=dict(color='darkgray', line_width=1, size=10),
                                 opacity=0.6, showlegend=False),
                      row=1, col=2)
        fig.add_trace(go.Scatter(x=df['Scaffold'], y=df['Buckling'], mode='markers', name='Buckling vs Estimated Scaffold Usage',
                                 marker=dict(color='darkgray', line_width=1, size=10), opacity=0.6,
                                 showlegend=False),
                      row=2, col=1)

        #### PLOTTING ALL PARETO FRONT:
        fig.add_trace(go.Scatter(x=pareto_rep_buck['Repulsion'], y=pareto_rep_buck['Buckling'], mode='markers',
                                 name='Buckling vs Repulsion', marker=dict(color='crimson', line_width=1, size=10),
                                 opacity=1.0, showlegend=False),
                      row=1, col=1)
        fig.add_trace(go.Scatter(x=pareto_rep_scaffold['Scaffold'], y=pareto_rep_scaffold['Repulsion'], mode='markers',
                                 name='Repulsion vs Estimated Scaffold Usage', marker=dict(color='crimson', line_width=1, size=10),
                                 opacity=1.0, showlegend=False),
                      row=1, col=2)
        fig.add_trace(go.Scatter(x=pareto_buck_scaffold['Scaffold'], y=pareto_buck_scaffold['Buckling'], mode='markers',
                                 name='Buckling vs Estimated Scaffold Usage',
                                 marker=dict(color='crimson', line_width=1, size=10), opacity=1.0,
                                 showlegend=False),
                      row=2, col=1)
        fig.add_trace(go.Scatter(x=[0.5], y=[0.5], text=['Epoch: 0 (post-archive initialization)'], mode="text",
                                 showlegend=False, textfont=dict(size=28)),
                      row=2, col=2)
        fig.update_xaxes(visible=False, row=2, col=2) # hide the axis and things
        fig.update_yaxes(visible=False, row=2, col=2)

        count = 1
        for i, j in zip([1, 1, 2], [1, 2, 1]):
            fig.update_xaxes(title_font=dict(size=20), tickfont=dict(size=16), row=i, col=j)
            fig.update_yaxes(title_font=dict(size=20), tickfont=dict(size=16), row=i, col=j)

            subplot_title = fig.layout.annotations[count - 1]['text']
            fig.update_annotations(font=dict(size=24), selector=dict(text=subplot_title))
            count += 1

        # Format of max_mins: [b_max, b_min, r_max, r_min, f_max, f_min]
        max_mins = find_min_and_maxes(l1=self.archive_post_temperature_initialization,
                                      l2=list(self.MOSA_archive.keys()))
        # Calculate padding so that points appear on the screen!
        padding_buckling = (max_mins[0] - max_mins[1]) * 0.1
        padding_repulsion = (max_mins[2] - max_mins[3]) * 0.1
        padding_scaffold = (max_mins[4] - max_mins[5]) * 0.1

        fig.update_xaxes(title_text='Repulsion', row=1, col=1, linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[3] - padding_repulsion, max_mins[2] + padding_repulsion])
        fig.update_yaxes(title_text='Buckling', row=1, col=1, linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[1] - padding_buckling, max_mins[0] + padding_buckling])

        fig.update_xaxes(title_text='Scaffold', row=1, col=2, linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[5] - padding_scaffold, max_mins[4] + padding_scaffold])
        fig.update_yaxes(title_text='Repulsion', row=1, col=2, linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[3] - padding_repulsion, max_mins[2] + padding_repulsion])

        fig.update_xaxes(title_text='Scaffold', row=2, col=1, linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[5] - padding_scaffold, max_mins[4] + padding_scaffold])
        fig.update_yaxes(title_text='Buckling', row=2, col=1, linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[1] - padding_buckling, max_mins[0] + padding_buckling])

        if save_or_show == 'show':
            fig.update_layout(title=dict(text='<b>' + plot_title + '</b>', font=dict(size=22, color='blue')))


        # Create Animation Frames:
        all_frames = []
        #all_frames_3D = []
        ct = 1
        fig.update_layout(plot_bgcolor='white')
        all_figs.append(fig)
        # Begin looping for all epochs
        for epoch, archive_at_epoch in self.pareto_animation.items():
            # Create the proper dataframe using the epoch data:
            df = pd.DataFrame(archive_at_epoch, columns=['Buckling', 'Repulsion', 'Scaffold'])
            df_repulsion_buckling = df[['Buckling', 'Repulsion']].copy()
            df_repulsion_scaffold = df[['Repulsion', 'Scaffold']].copy()
            df_buckling_scaffold = df[['Buckling', 'Scaffold']].copy()
            #df_full = df.copy()

            # Find the pareto optimal points:
            pareto_rep_buck = self.calculate_pareto_front(df=df_repulsion_buckling, field1='Repulsion',
                                                          field2='Buckling')
            pareto_rep_scaffold = self.calculate_pareto_front(df=df_repulsion_scaffold, field1='Repulsion', field2='Scaffold')
            pareto_buck_scaffold = self.calculate_pareto_front(df=df_buckling_scaffold, field1='Buckling', field2='Scaffold')
            #pareto_3D = self.calculate_3D_pareto(df=df_full)

            # Create the Data traces we are updating:
            data = [
                go.Scatter(x=df['Repulsion'], y=df['Buckling'], mode='markers',
                           name='Buckling vs Repulsion', marker=dict(color='darkgray', line_width=1, size=10),
                           opacity=0.6, showlegend=False),
                go.Scatter(x=df['Scaffold'], y=df['Repulsion'], mode='markers', name='Repulsion vs Estimated Scaffold Usage',
                           marker=dict(color='darkgray', line_width=1, size=10),
                           opacity=0.6, showlegend=False),
                go.Scatter(x=df['Scaffold'], y=df['Buckling'], mode='markers', name='Buckling vs Estimated Scaffold Usage',
                           marker=dict(color='darkgray', line_width=1, size=10), opacity=0.6,
                           showlegend=False),
                go.Scatter(x=pareto_rep_buck['Repulsion'], y=pareto_rep_buck['Buckling'], mode='markers',
                          name='Buckling vs Repulsion', marker=dict(color='crimson', line_width=1, size=10),
                          opacity=1.0, showlegend=False),
                go.Scatter(x=pareto_rep_scaffold['Scaffold'], y=pareto_rep_scaffold['Repulsion'], mode='markers',
                           name='Repulsion vs Estimated Scaffold Usage', marker=dict(color='crimson', line_width=1, size=10),
                           opacity=1.0, showlegend=False),
                go.Scatter(x=pareto_buck_scaffold['Scaffold'], y=pareto_buck_scaffold['Buckling'], mode='markers',
                           name='Buckling vs Estimated Scaffold Usage',
                           marker=dict(color='crimson', line_width=1, size=10), opacity=1.0,
                           showlegend=False),
                go.Scatter(x=[0.5], y=[0.5], text=[f'Epoch: {ct} (post-archive initialization)'], mode="text",
                           textfont=dict(size=28), showlegend=False)
            ]

            if save_or_show == 'show':
                traces = [0, 1, 2, 3, 4, 5, 6]
                #traces_2 = [0, 1]
                # Append to all_frames:
                all_frames.append(dict(name=epoch, data=data, traces=traces))
                #all_frames_3D.append(dict(name=epoch, data=data2, traces=traces_2))
            else:
                new_fig = make_subplots(rows=2, cols=2,
                            subplot_titles=("Repulsion vs Buckling", "Repulsion vs Estimated Scaffold Usage",
                                            "Buckling vs Estimated Scaffold Usage", ""),
                            start_cell="top-left",
                            specs=[
                                [{"type": "scatter"}, {"type": "scatter"}],
                                [{"type": "scatter"}, {"type": "scatter"}]
                            ]
                            )

                new_fig.add_trace(data[0], row=1, col=1)
                new_fig.add_trace(data[1], row=1, col=2)
                new_fig.add_trace(data[2], row=2, col=1)
                new_fig.add_trace(data[3], row=1, col=1)
                new_fig.add_trace(data[4], row=1, col=2)
                new_fig.add_trace(data[5], row=2, col=1)
                new_fig.add_trace(data[6], row=2, col=2)

                # Update layout of new_fig
                new_fig.update_xaxes(visible=False, row=2, col=2)  # hide the axis and things
                new_fig.update_yaxes(visible=False, row=2, col=2)
                count = 1
                for i, j in zip([1, 1, 2], [1, 2, 1]):
                    new_fig.update_xaxes(title_font=dict(size=20), tickfont=dict(size=16), row=i, col=j)
                    new_fig.update_yaxes(title_font=dict(size=20), tickfont=dict(size=16), row=i, col=j)

                    subplot_title = new_fig.layout.annotations[count - 1]['text']
                    new_fig.update_annotations(font=dict(size=24), selector=dict(text=subplot_title))
                    count += 1

                # Format of max_mins: [b_max, b_min, r_max, r_min, f_max, f_min]
                max_mins = find_min_and_maxes(l1=self.archive_post_temperature_initialization,
                                              l2=list(self.MOSA_archive.keys()))
                # Calculate padding so that points appear on the screen!
                padding_buckling = (max_mins[0] - max_mins[1]) * 0.1
                padding_repulsion = (max_mins[2] - max_mins[3]) * 0.1
                padding_scaffold = (max_mins[4] - max_mins[5]) * 0.1

                new_fig.update_xaxes(title_text='Repulsion', row=1, col=1, linecolor='black', linewidth=1, mirror=True,)
                                 #range=[max_mins[3] - padding_repulsion, max_mins[2] + padding_repulsion])
                new_fig.update_yaxes(title_text='Buckling', row=1, col=1, linecolor='black', linewidth=1, mirror=True,)
                                 #range=[max_mins[1] - padding_buckling, max_mins[0] + padding_buckling])

                new_fig.update_xaxes(title_text='Estimated Scaffold Usage', row=1, col=2, linecolor='black', linewidth=1, mirror=True,)
                                 #range=[max_mins[5] - padding_scaffold, max_mins[4] + padding_scaffold])
                new_fig.update_yaxes(title_text='Repulsion', row=1, col=2, linecolor='black', linewidth=1, mirror=True,)
                                 #range=[max_mins[3] - padding_repulsion, max_mins[2] + padding_repulsion])

                new_fig.update_xaxes(title_text='Estimated Scaffold Usage', row=2, col=1, linecolor='black', linewidth=1, mirror=True,)
                                 #range=[max_mins[5] - padding_scaffold, max_mins[4] + padding_scaffold])
                new_fig.update_yaxes(title_text='Buckling', row=2, col=1, linecolor='black', linewidth=1, mirror=True,)
                                 #range=[max_mins[1] - padding_buckling, max_mins[0] + padding_buckling])

                new_fig.update_layout(plot_bgcolor='white')
                all_figs.append(new_fig)
                traces = [0, 1, 2, 3, 4, 5, 6]
                # traces_2 = [0, 1]
                # Append to all_frames:
                all_frames.append(dict(name=epoch, data=data, traces=traces))
                del new_fig

            del df  # Clear out df for memory purposes...
            ct += 1

        # Creating menus and slider bars in plotly:
        updatemenus = [dict(type='buttons',
                            buttons=[dict(label='Play',
                                          method='animate',
                                          args=[None,
                                                dict(frame=dict(duration=600, redraw=True),
                                                     transition=dict(duration=0),
                                                     fromcurrent=True,
                                                     mode='immediate'
                                                     )]),
                                     dict(label='Pause',
                                          method='animate',
                                          args=[[None],
                                                dict(frame=dict(duration=0, redraw=False),
                                                     mode='immediate',
                                                     transition=dict(duration=0)
                                                     )]
                                          )
                                     ],
                            direction='left',
                            pad=dict(r=10, t=85),
                            showactive=True, x=0.1, y=0, xanchor='right', yanchor='top')
                       ]
        fig.update(frames=all_frames)
        fig.update_layout(plot_bgcolor='white')
        fig.update_layout(updatemenus=updatemenus)
        if save_or_show == 'show':
            # sliders=sliders)
            # fig_3D.update(frames=all_frames_3D)
            # fig_3D.update_layout(plot_bgcolor='white')
            # fig_3D.update_layout(updatemenus=updatemenus)
            #self.all_figures_for_animation = all_figs
            fig.show()
        else:
            self.all_figures_for_animation = [fig]  # Save this for future loading purposes
            return all_figs


    def explore_designs_from_final_archive(self, archive_data=None) -> list:
        """
        This function will update a 3D scatter plot showing the designs that were generated during a MOSA simulation.
        Currently, they will just be a slider on an updating plotly figure. It would be cool to add an "export PLY"
        button that automatically exports the PLY file (or DNA files) for the currently selected figure...

        The user could pass in some archive_data dictionary to use, otherwise the final archive stored in MOSA_archive
        will be exported for this plot.
        """
        if archive_data is None:
            archive_data = self.MOSA_archive

        figure_list = []
        for objectives, design in archive_data.items():
            figure_list.append(design.display_graph(return_figure=True))

        # Return the Figure List:
        return figure_list


    def paretos_for_dash_app(self, data_points: list, max_mins: list):
        """
        This function will use plotly to plot (and save) the pareto front. We save out the Pareto at the end of the
        initial temperature optimization during White's algorithm and at the end of the MOSA process.
        """
        # We will create side-by-side plots so we can see how the pareto optimal surface came to be during MOSA:
        rep_vs_buck = go.Figure()
        rep_vs_scaffold = go.Figure()
        buck_vs_scaffold = go.Figure()

        # Process data to add a singular trace to the figure to make life easier:
        df = pd.DataFrame(data_points, columns=['Buckling', 'Repulsion', 'Scaffold'])
        df_repulsion_buckling = df[['Repulsion', 'Buckling']].copy()
        df_repulsion_scaffold = df[['Repulsion', 'Scaffold']].copy()
        df_buckling_scaffold = df[['Buckling', 'Scaffold']].copy()

        # Find the pareto optimal points:
        pareto_rep_buck = self.calculate_pareto_front(df=df_repulsion_buckling, field1='Buckling', field2='Repulsion')
        pareto_rep_scaffold = self.calculate_pareto_front(df=df_repulsion_scaffold, field1='Repulsion', field2='Scaffold')
        pareto_buck_scaffold = self.calculate_pareto_front(df=df_buckling_scaffold, field1='Buckling', field2='Scaffold')

        #### PLOTTING ALL PARETO POINTS (NOT THE PARETO FRONT):
        rep_vs_buck.add_trace(go.Scatter(x=df['Buckling'], y=df['Repulsion'], mode='markers',
                                 name='Repulsion vs Buckling', marker=dict(color='darkgray', line_width=1, size=10),
                                 opacity=0.6, showlegend=False))
        rep_vs_scaffold.add_trace(go.Scatter(x=df['Scaffold'], y=df['Repulsion'], mode='markers', name='Repulsion vs Estimated Scaffold Usage',
                                 marker=dict(color='darkgray', line_width=1, size=10),
                                 opacity=0.6, showlegend=False))
        buck_vs_scaffold.add_trace(go.Scatter(x=df['Scaffold'], y=df['Buckling'], mode='markers', name='Buckling vs Estimated Scaffold Usage',
                                 marker=dict(color='darkgray', line_width=1, size=10), opacity=0.6,
                                 showlegend=False))

        #### PLOTTING ALL PARETO FRONT:
        rep_vs_buck.add_trace(go.Scatter(x=pareto_rep_buck['Buckling'], y=pareto_rep_buck['Repulsion'], mode='markers',
                                 name='Repulsion vs Buckling', marker=dict(color='crimson', line_width=1, size=10),
                                 opacity=1.0, showlegend=False))
        rep_vs_scaffold.add_trace(go.Scatter(x=pareto_rep_scaffold['Scaffold'], y=pareto_rep_scaffold['Repulsion'], mode='markers',
                                 name='Repulsion vs Estimated Scaffold Usage', marker=dict(color='crimson', line_width=1, size=10),
                                 opacity=1.0, showlegend=False))

        buck_vs_scaffold.add_trace(go.Scatter(x=pareto_buck_scaffold['Scaffold'], y=pareto_buck_scaffold['Buckling'], mode='markers',
                                 name='Buckling vs Estimated Scaffold Usage',
                                 marker=dict(color='crimson', line_width=1, size=10), opacity=1.0,
                                 showlegend=False))

        for fig in [rep_vs_buck, rep_vs_scaffold, buck_vs_scaffold]:
            fig.update_xaxes(title_font=dict(size=20), tickfont=dict(size=16))
            fig.update_yaxes(title_font=dict(size=20), tickfont=dict(size=16))


        # Format of max_mins: [b_max, b_min, r_max, r_min, f_max, f_min]
        # Calculate padding so that points appear on the screen!
        padding_buckling = (max_mins[0] - max_mins[1]) * 0.1
        padding_repulsion = (max_mins[2] - max_mins[3]) * 0.1
        padding_scaffold = (max_mins[4] - max_mins[5]) * 0.1

        rep_vs_buck.update_xaxes(title_text='Buckling', linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[1] - padding_buckling, max_mins[0] + padding_buckling])
        rep_vs_buck.update_yaxes(title_text='Repulsion', linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[3] - padding_repulsion, max_mins[2] + padding_repulsion])

        rep_vs_scaffold.update_xaxes(title_text='Estimated Scaffold Usage', linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[5] - padding_scaffold, max_mins[4] + padding_scaffold])
        rep_vs_scaffold.update_yaxes(title_text='Repulsion', linecolor='black', linewidth=1, mirror=True,)
                         #range=[max_mins[3] - padding_repulsion, max_mins[2] + padding_repulsion])

        buck_vs_scaffold.update_xaxes(title_text='Estimated Scaffold Usage', linecolor='black', linewidth=1, mirror=True)
                         #range=[max_mins[5] - padding_scaffold, max_mins[4] + padding_scaffold])
        buck_vs_scaffold.update_yaxes(title_text='Buckling', linecolor='black', linewidth=1, mirror=True)
                         #range=[max_mins[1] - padding_buckling, max_mins[0] + padding_buckling])

        rep_vs_buck.update_layout(plot_bgcolor='white')
        rep_vs_scaffold.update_layout(plot_bgcolor='white')
        buck_vs_scaffold.update_layout(plot_bgcolor='white')

        return rep_vs_buck, rep_vs_scaffold, buck_vs_scaffold


    """def update_pareto_color(self, data_points: list, max_mins: list, f1_x: str, f2_y: str, changePt: tuple):
        # We will create side-by-side plots so we can see how the pareto optimal surface came to be during MOSA:
        updated_fig = go.Figure()

        # Process data to add a singular trace to the figure to make life easier:
        df = pd.DataFrame(data_points, columns=['Buckling', 'Repulsion', 'Faces'])
        df_desired = df[[f1_x, f2_y]].copy()
        mask = df.apply(lambda row: tuple(row) != changePt, axis=1)
        df = df[mask]

        # Find the pareto optimal points:
        pareto_points = self.calculate_pareto_front(df=df_desired, field1=f1_x, field2=f2_y)
        mask = pareto_points.apply(lambda row: tuple(row) != changePt, axis=1)
        pareto_points = pareto_points[mask]

        #### PLOTTING ALL PARETO POINTS (NOT THE PARETO FRONT):
        updated_fig.add_trace(go.Scatter(x=df[f1_x], y=df[f2_y], mode='markers',
                                         name='Repulsion vs Buckling',
                                         marker=dict(color='darkgray', line_width=1, size=10),
                                         opacity=0.6, showlegend=False))

        #### PLOTTING ALL PARETO FRONT:
        updated_fig.add_trace(go.Scatter(x=pareto_points[f1_x], y=pareto_points[f2_y], mode='markers',
                                         name='Repulsion vs Buckling',
                                         marker=dict(color='crimson', line_width=1, size=10),
                                         opacity=1.0, showlegend=False))

        ### PLOTTING THE ACTIVE POINT:
        updated_fig.add_trace(go.Scatter(x=[changePt[0]], y=[changePt[1]], mode='markers',
                                         marker=dict(color='green', line_width=1, size=10, symbol='square'),
                                         opacity=1.0, showlegend=False))


        updated_fig.update_xaxes(title_font=dict(size=20), tickfont=dict(size=16))
        updated_fig.update_yaxes(title_font=dict(size=20), tickfont=dict(size=16))

        # Format of max_mins: [b_max, b_min, r_max, r_min, f_max, f_min]
        # Calculate padding so that points appear on the screen!
        padding_buckling = (max_mins[0] - max_mins[1]) * 0.1
        padding_repulsion = (max_mins[2] - max_mins[3]) * 0.1
        padding_faces = (max_mins[4] - max_mins[5]) * 0.1

        if f1_x == 'Buckling' and f2_y == 'Repulsion':
            updated_fig.update_xaxes(title_text='Buckling', linecolor='black', linewidth=1, mirror=True,
                                     range=[max_mins[1] - padding_buckling, max_mins[0] + padding_buckling])
            updated_fig.update_yaxes(title_text='Repulsion', linecolor='black', linewidth=1, mirror=True,
                                     range=[max_mins[3] - padding_repulsion, max_mins[2] + padding_repulsion])
        elif f1_x == 'Faces' and f2_y == 'Repulsion':
            updated_fig.update_xaxes(title_text='Faces', linecolor='black', linewidth=1, mirror=True,
                                      range=[max_mins[5] - padding_faces, max_mins[4] + padding_faces])
            updated_fig.update_yaxes(title_text='Repulsion', linecolor='black', linewidth=1, mirror=True,
                                      range=[max_mins[3] - padding_repulsion, max_mins[2] + padding_repulsion])
        else:
            updated_fig.update_xaxes(title_text='Faces', linecolor='black', linewidth=1, mirror=True,
                                       range=[max_mins[5] - padding_faces, max_mins[4] + padding_faces])
            updated_fig.update_yaxes(title_text='Buckling', linecolor='black', linewidth=1, mirror=True,
                                       range=[max_mins[1] - padding_buckling, max_mins[0] + padding_buckling])

        updated_fig.update_layout(plot_bgcolor='white')

        return updated_fig"""

