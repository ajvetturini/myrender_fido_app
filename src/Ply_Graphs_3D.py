"""
AJ Vetturini
IDIG and MMBL
Carnegie Mellon University
Advised By: Jon Cagan and Rebecca Taylor

This program is used to support the Shape Annealing of 3D PLY Files which will be used for DNA Origami Generation along
with Shape_Annealing_3DPly, and built into one large package which will be released.
"""

# Import Modules:
from dataclasses import dataclass, field
from collections import defaultdict, Counter
from typing import Optional, Union
import numpy as np
import networkx as nx
import plotly.graph_objs as go
import os
from scipy.spatial import Delaunay, ConvexHull
from sklearn.decomposition import PCA
import itertools
from random import choice, seed, randint
import warnings
from copy import deepcopy
import trimesh
import pickle
import plotly.offline as pyo
import plotly.io as pio
import json

# Global Variables:
numDecimals = 3  # This is used when rounding throughout the program to calculate the number of decimals used in calcs.
title_font_size = 24
axis_font_size = 32


# Dataclass for the Nanoparticles / areas to "avoid" when creating a 3D Design Space. Since the edges matter for
# calculating amount of DNA used, I am initially only using 'box' zones. Might add a new type later.
@dataclass
class NPBox(object):
    # We locate the Nanoparticle based on the bottom left corner and the top right corners 3D Coordinates as a 1D array
    # These MUST be OPPOSING corners to define out the "bounding box" of sorts.
    c1: np.array
    c2: np.array
    relabel_dict: dict = field(default_factory=dict)
    faces: dict = field(default_factory=dict)

    def __post_init__(self):
        # First we call the verify function to ensure that these point ARE in fact opposing:
        self.verify_corners_opposing()

        # Create a tri-mesh box:
        self.cube_mesh = trimesh.creation.box(bounds=np.array([self.c1, self.c2]))

        self.NP_as_graph = nx.Graph()  # Initialize a graph that will store data of these points
        # Calculate the differences in X Y and Z:
        global numDecimals
        dx = round(self.c2[0] - self.c1[0], numDecimals)
        dy = round(self.c2[1] - self.c1[1], numDecimals)
        dz = round(self.c2[2] - self.c1[2], numDecimals)

        # Now with these values defined, we can find the other coordinates for the graph as follows:
        graph_points = [(0, {'x': self.c1[0], 'y': self.c1[1], 'z': self.c1[2], 'terminal': True}),
                        (1, {'x': self.c1[0] + dx, 'y': self.c1[1], 'z': self.c1[2], 'terminal': True}),
                        (2, {'x': self.c1[0], 'y': self.c1[1] + dy, 'z': self.c1[2], 'terminal': True}),
                        (3, {'x': self.c1[0], 'y': self.c1[1], 'z': self.c1[2] + dz, 'terminal': True}),
                        (4, {'x': self.c2[0] - dx, 'y': self.c2[1], 'z': self.c2[2], 'terminal': True}),
                        (5, {'x': self.c2[0], 'y': self.c2[1] - dy, 'z': self.c2[2], 'terminal': True}),
                        (6, {'x': self.c2[0], 'y': self.c2[1], 'z': self.c2[2] - dz, 'terminal': True}),
                        (7, {'x': self.c2[0], 'y': self.c2[1], 'z': self.c2[2], 'terminal': True})
                        ]
        graph_edges = [(0, 1), (0, 2), (0, 3), (3, 4), (3, 5), (4, 2), (1, 5), (1, 6), (6, 2), (6, 7), (7, 5), (7, 4)]

        # Creating the graph from these:
        self.NP_as_graph.add_nodes_from(graph_points)
        self.NP_as_graph.add_edges_from(graph_edges)

        # Create faces which will be read in via the DesignSpace:
        all_faces = [(2, 6, 1, 0), (2, 0, 3, 4), (4, 3, 5, 7), (7, 5, 1, 6), (2, 4, 7, 6), (0, 1, 5, 3)]
        for f in all_faces:
            # We can't determine face-directionality until the 3D structure is finalized. However, there are certain
            # functionalities we will want for the face:
            verts = []
            for v in list(f):
                # Grab the X, Y, and Z values:
                pt = [self.NP_as_graph.nodes[v]['x'], self.NP_as_graph.nodes[v]['y'], self.NP_as_graph.nodes[v]['z']]
                verts.append(np.array(pt))
            self.faces[f] = MeshFaces(v1=verts[0], v2=verts[1], v3=verts[2], v4=verts[3], face_type='rect')

        # Append a list:
        self.graph_points = []
        for i in graph_points:
            self.graph_points.append(i[1])

    def point_inside_space(self, pt: np.array) -> bool:
        """
        This function will determine if a passed in point is located inside this space, and therefore an invalid
        geometry!!
        :param pt: XYZ of a point we are checking
        :return: True: Point inside of the NP space; False: point not inside
        """
        if (self.c1[0] < pt[0] < self.c2[0] or self.c2[0] < pt[0] < self.c1[0]) and \
                (self.c1[1] < pt[1] < self.c2[1] or self.c2[1] < pt[1] < self.c1[1]) and \
                (self.c1[2] < pt[2] < self.c2[2] or self.c2[2] < pt[2] < self.c1[2]):
            return True
        # If we are not within the box, we return False
        return False

    def edge_intersecting_face(self, p1: np.array, p2: np.array) -> bool:
        """
        This function will determine if a passed in line segment intersects a face of an NP Region.
        :param p1: Point 1 of a line segment definition
        :param p2: Point 2 of a line segment definition
        :return: True: Edge intersects the NP space; False: edge does not intersect
        """
        face_intersections = []
        u = p2 - p1
        intersector = trimesh.ray.ray_triangle.RayMeshIntersector(self.cube_mesh)
        intersections, ignore1, ignore2 = intersector.intersects_location([p1], [u], multiple_hits=False)

        # Now we check to verify the intersections:
        for hit in intersections:
            if np.array_equal(p1, hit) or np.array_equal(p2, hit):
                pass
            else:
                # In this case, we also need to check to see if this intersection hit is inbetween our two points or
                # not, since a Ray does not have any bounds.
                d_12 = np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2)
                d_1_hit = np.sqrt((hit[0] - p1[0])**2 + (hit[1] - p1[1])**2 + (hit[2] - p1[2])**2)
                d_2_hit = np.sqrt((p2[0] - hit[0])**2 + (p2[1] - hit[1])**2 + (p2[2] - hit[2])**2)
                # Now, if the following check is true, then we have a face intersection:
                if d_1_hit < d_12 and d_2_hit < d_12:
                    face_intersections.append(True)

        if any(face_intersections):
            return True
        else:
            return False

    def verify_corners_opposing(self) -> None:
        """
        This function verifies that the corners passed in ARE in fact opposing corners
        :return: Nothing, raises an exception if they are invalid.
        """
        # Basically, all values of C1 MUST be less than C2:
        if any(C1 >= C2 for C1, C2 in zip(self.c1, self.c2)):
            raise ValueError('The input conditions for the opposing corners must be that all values in c1 are less than'
                             ' the values found in c2 for all of X, Y, Z.')


# These are faces that will be stored during the 3D Space Generation
@dataclass
class MeshFaces(object):
    # Our faces can be either a triangle or a square, but mainly a triangle.
    v1: np.array
    v2: np.array
    v3: np.array
    v4: Optional[np.array] = None
    faces: Optional[list] = None
    face_type: str = 'triangle'
    area: float = None


    def __post_init__(self):
        if self.face_type != 'triangle' and self.face_type != 'rect':
            raise Exception(f'Invalid face type passed in: {self.face_type}')
        self.face_angles = [90, 90, 90, 90]
        if self.face_type == 'triangle':
            #self.calculate_face_angles()
            # Calculate the area:
            AB = self.v2 - self.v1
            AC = self.v3 - self.v1
            self.area = 0.5 * np.linalg.norm(np.cross(AB, AC))
        elif self.face_type == 'rect':
            # Calculating area of a rectangular face by adding two triangles:
            AB = self.v2 - self.v1
            AC = self.v3 - self.v1
            A1 = 0.5 * np.linalg.norm(np.cross(AB, AC))

            AD = self.v4 - self.v1
            A2 = 0.5 * np.linalg.norm(np.cross(AC, AD))
            self.area = A1 + A2

            pass


    def calculate_face_angles(self) -> None:
        """
        This function is called during post_init to calculate the face angles for the vertices. This is only called
        for the triangular faces!
        :return: Nothing
        """
        AB = self.v2 - self.v1
        BC = self.v3 - self.v2
        CA = self.v1 - self.v3
        # Calculate Magnitudes and Products:
        mag_AB = np.linalg.norm(AB)
        mag_BC = np.linalg.norm(BC)
        mag_CA = np.linalg.norm(CA)

        # Calculate angles using dot product formula
        angle_ABC = np.arccos(np.dot(AB, -BC) / (mag_AB * mag_BC))
        angle_BCA = np.arccos(np.dot(BC, -CA) / (mag_BC * mag_CA))
        angle_CAB = np.arccos(np.dot(CA, -AB) / (mag_CA * mag_AB))

        # Convert angles to degrees
        global numDecimals
        angle_ABC_deg = np.round(np.degrees(angle_ABC), numDecimals)
        angle_BCA_deg = np.round(np.degrees(angle_BCA), numDecimals)
        angle_CAB_deg = np.round(np.degrees(angle_CAB), numDecimals)

        # Reassign the face angles:
        self.face_angles = [angle_ABC_deg, angle_BCA_deg, angle_CAB_deg]

    def review_vertex_orientation(self, reference_point: np.array):
        """
        This function is use to hopefully try and ensure that these MeshFaces dataclass are always CCW orientation
        :param reference_point: A point passed in to serve as a reference to establish "directionality" for calculating
                                plane values (a, b, c, d)
        :return: Calculated values of a, b, c, and d
        """
        ## FIRST: Calculate face normal with the current orientation:
        if self.face_type == 'triangle':
            vertices = [self.v1, self.v2, self.v3]
        else:
            vertices = [self.v1, self.v2, self.v3, self.v4]

        global numDecimals
        edge1 = np.round((vertices[1] - vertices[0]), numDecimals)
        edge2 = np.round((vertices[-1] - vertices[0]), numDecimals)
        face_normal = np.cross(edge1, edge2)

        # THIRD: Calculate thee centroid and dot product:
        centroid_to_reference = reference_point - np.mean(vertices, axis=0)
        if np.dot(face_normal, centroid_to_reference) < 0:
            # In this case, the vertices are NOT CCW and MUST be re-arranged before calculating ABCD
            if self.face_type == 'triangle':
                correct_order = [self.v1, self.v3, self.v2]
            else:
                correct_order = [self.v1, self.v4, self.v3, self.v2]
        else:
            # Otherwise, the values ARE in the correct order so no need to change:
            if self.face_type == 'triangle':
                correct_order = [self.v1, self.v2, self.v3]
            else:
                correct_order = [self.v1, self.v2, self.v3, self.v4]

        # FOURTH AND FINALLY: Calculate A, B, C, D values for the plane:
        AB = correct_order[1] - correct_order[0]
        AC = correct_order[-1] - correct_order[0]
        norm = np.cross(AB, AC)
        nD = 0
        return np.round(norm, nD)

    def equation_of_face(self) -> tuple[float, float, float, float]:
        """
        This function will calculate the equation of the triangular face. This function is solely called by
        find_coplanar_faces. This is why I don't care what the exact orientation of the face normals is because I only
        care for 180 degrees between planes so it doesn't matter where the angle is taken from.

        :return: Values of a, b, c, and d to define a plane
        """
        AB = self.v2 - self.v1
        AC = self.v3 - self.v1
        a, b, c = np.cross(AB, AC)
        d = -1 * np.dot(np.cross(AB, AC), self.v1)
        return a, b, c, d

    def get_face_plot(self, color_type: str, override_opacity: bool = False) -> go.Mesh3d:
        """
        This function will be used to plot the mesh face simply via Plotly
        :param color_type: Dictates which color to use in the plot via the plotting dictionary
        :return: A Mesh3D value that will be used to plot
        """
        # Color dictionary used in plotting the faces below:
        color_key = {'NP': 'red', 'binding_region': 'green', 'other': 'grey'}

        if self.face_type == 'rect':
            all_verts = np.stack((self.v1, self.v2, self.v3, self.v4), axis=0)
            i = [0, 0]
            j = [1, 2]
            k = [2, 3]
        else:
            all_verts = np.stack((self.v1, self.v2, self.v3), axis=0)
            i = [0]
            j = [1]
            k = [2]
        # Now create the plot trace and return:
        if override_opacity:
            return go.Mesh3d(x=all_verts[:, 0], y=all_verts[:, 1], z=all_verts[:, 2], i=i, j=j, k=k,
                             color=color_key[color_type], flatshading=True, opacity=0.9)
        else:
            return go.Mesh3d(x=all_verts[:, 0], y=all_verts[:, 1], z=all_verts[:, 2], i=i, j=j, k=k,
                             color=color_key[color_type], flatshading=True, opacity=1)

    def divide_face(self) -> tuple[np.array, np.array, tuple, list, list]:
        """
        This function will divide a mesh face into two by selecting a node and bisecting it to the far edge.
        :return: It will return 1) Tuple(newEdgePoint, NewPointPosition, DividedEdgeVerts, NewFace1, NewFace2)
        """
        # We want to divide the triangular face in such a way that gurantees the edge we create is
        all_verts = [self.v1, self.v2, self.v3]
        randIndex = randint(0, 2)  # Select a random vertex
        A = all_verts[randIndex]
        # Next we calculate the midpoints from A to the opposite index
        B = (A + all_verts[(randIndex + 1) % 3]) / 2

        # Finally, we define the new triangles based on XYZ Positions:
        newTriangle1 = [A, B, all_verts[(randIndex + 2) % 3]]
        newTriangle2 = [all_verts[(randIndex + 1) % 3], B, all_verts[(randIndex + 2) % 3]]

        # Now return these points to create / add new faces:
        return all_verts[(randIndex + 2) % 3], B, (A, all_verts[(randIndex + 1) % 3]), newTriangle1, newTriangle2

    def update_vertex_values(self, old_xyz: np.array, new_xyz: np.array) -> None:
        """
        This function is called whenever the vertices of a face change values. This function will also update teh abcd
        values.
        :param old_xyz: Old position to search for to change the correct vertex XYZ
        :param new_xyz: New position to update with
        :return: Nothing, just changing values of object
        """
        # We first ensure decimal place values:
        global numDecimals
        rounded_new_xyz = np.round(new_xyz, decimals=numDecimals)

        # I do not know if there is a better way to do this, but since the faces only have a few values I am just
        # going to use if statements. I tried a for loop but that didn't update the object value:
        if np.array_equal(old_xyz, self.v1):
            self.v1 = rounded_new_xyz
        elif np.array_equal(old_xyz, self.v2):
            self.v2 = rounded_new_xyz
        elif np.array_equal(old_xyz, self.v3):
            self.v3 = rounded_new_xyz
        elif np.array_equal(old_xyz, self.v4):
            self.v4 = rounded_new_xyz
        else:
            raise Exception('Unable to find original vertex value and therefore there is a logic error --> ERROR!!!')

    def calculate_face_centroid(self) -> np.array:
        """
        This function calculates and returns the centroid of the face to use specifically with the Retriangulate Rule
        :return: XYZ Coordinates of Face Centroid
        """
        global numDecimals
        if self.face_type == 'triangle':
            return np.round(((self.v1 + self.v2 + self.v3) / 3), decimals=numDecimals)
        else:
            return np.round(((self.v1 + self.v2 + self.v3 + self.v4) / 4), decimals=numDecimals)

    def return_both_normal_vectors(self, reference_pt: np.array) -> np.array:
        ## FIRST: Calculate face normal with the current orientation:
        if self.face_type == 'triangle':
            vertices = [self.v1, self.v2, self.v3]
        else:
            vertices = [self.v1, self.v2, self.v3, self.v4]

        global numDecimals
        edge1 = np.round((vertices[1] - vertices[0]), numDecimals)
        edge2 = np.round((vertices[-1] - vertices[0]), numDecimals)
        face_normal = np.cross(edge1, edge2)

        # THIRD: Calculate thee centroid and dot product:
        centroid_to_reference = reference_pt - np.mean(vertices, axis=0)
        if np.dot(face_normal, centroid_to_reference) < 0:
            # In this case, the vertices are NOT CCW and MUST be re-arranged before calculating ABCD
            if self.face_type == 'triangle':
                correct_order = [self.v1, self.v3, self.v2]
            else:
                correct_order = [self.v1, self.v4, self.v3, self.v2]
        else:
            # Otherwise, the values ARE in the correct order so no need to change:
            if self.face_type == 'triangle':
                correct_order = [self.v1, self.v2, self.v3]
            else:
                correct_order = [self.v1, self.v2, self.v3, self.v4]

        # Once the correct order is established, we can quickly find the normal vectors:
        edge1 = correct_order[1] - correct_order[0]
        edge2 = correct_order[2] - correct_order[0]
        normal = np.cross(edge1, edge2)
        normal = normal / np.linalg.norm(normal)
        return normal, -1 * normal

    def recalculate_face_area(self) -> float:
        AB = self.v2 - self.v1
        AC = self.v3 - self.v1
        self.area = 0.5 * np.linalg.norm(np.cross(AB, AC))
        return self.area


# Dataclass that will represent either a binding vertex, a binding edge, or a binding face for our 3D Ply Generation
@dataclass
class BindingRegion(object):
    # A binding vertex will only have one passed in point, an edge will have two points, and a face will have 3 points.
    # A binding face can also have a 4th index for a rectangular face. This will create a MeshFaces objective if 3 or 4.
    v1: np.array
    v2: Optional[np.array] = None
    v3: Optional[np.array] = None
    v4: Optional[np.array] = None
    relabel_dict: dict = field(default_factory=dict)

    def __post_init__(self):
        # First we check to see what kind of face we have:
        if self.v2 is None and self.v3 is None and self.v4 is None:
            self.binding_type = 'vertex'
            v_list = [self.v1]
            graph_edges = []
        elif self.v3 is None and self.v4 is None:
            self.binding_type = 'edge'
            v_list = [self.v1, self.v2]
            graph_edges = [(0, 1)]
        elif self.v4 is None:
            self.binding_type = 'triangle'
            v_list = [self.v1, self.v2, self.v3]
            graph_edges = [(0, 1), (0, 2), (1, 2)]
        elif self.v4 is not None:
            self.binding_type = 'rect'
            v_list = [self.v1, self.v2, self.v3, self.v4]
            graph_edges = [(0, 1), (0, 3), (1, 2), (2, 3)]
        else:
            raise Exception('Invalid entries for the binding region')

        # Next we create a NetworkX Graph:
        self.binder_as_graph = nx.Graph()  # Initialize a graph that will store data of these points
        graph_points = []
        for ct, i in enumerate(v_list):
            graph_points.append((ct, {'x': i[0], 'y': i[1], 'z': i[2], 'terminal': True, 'binding_zone': True}))

        # Now add nodes and points to the NetworkX Graph:
        self.binder_as_graph.add_nodes_from(graph_points)
        self.binder_as_graph.add_edges_from(graph_edges)

        # Append a list for all the points found in the binding_region:
        self.graph_points = []
        for i in graph_points:
            self.graph_points.append(i[1])

        # Now if we are a triangle or rect binding type we create a face:
        if self.binding_type == 'triangle':
            self.face = ([0, 1, 2], MeshFaces(v1=self.v1, v2=self.v2, v3=self.v3, face_type='triangle'))
        elif self.binding_type == 'rect':
            self.face = ([0, 1, 2, 3], MeshFaces(v1=self.v1, v2=self.v2, v3=self.v3, v4=self.v4, face_type='rect'))
        else:
            self.face = None


## GLOBAL VARIABLES:
# This is a variable used in design constraint checking to see if an applied rule is invalid
design_constraint_failure = False

# Master Dataclass that will represent the changing design space and utilize the above as inputs.
@dataclass
class DesignSpace3D(object):
    # User Defined Variables (must be passed in)
    maxX: float
    maxY: float
    maxZ: float
    NP_Regions: list
    binding_regions: list
    min_edge_length: int  # This is in units of BASEPAIRS and can NOT be less than 31!

    # Preset Variables that users can modify if desired / they know what they are doing:
    max_edge_length_multiplier: float = 5  # This is a value of N*min_edge_length for a max edge length. This value is N
    random_seed: int = 520
    routing_software: str = 'DAEDALUS'
    extend_rule_distance: float = 1.5  # In units of nm. This can be easily changed, 1 seems to work fine.
    max_number_basepairs: int = 7249  # Default value for M13 phage
    min_face_angle: float = 8.  # This is the minimal angle a face can have IN DEGREES!!!
    persistance_length: float = 45.  # This is the persistance length of DNA by default but can be changed by user.
    dna_diameter: float = 2.  # Diameter of DNA in units of nm
    max_number_edges_at_node: int = 6
    check_scaffold_length_constraint: bool = True  # Dictates if the structures we generate "care" about scaffold length
    allow_self_intersecting_geometry: bool = True  # Dictates if we allow invalid mesh geometry to appear or not
    repulsion_distance: float = 2.  # This is the value of distance between centers of edges to "count" repulsion forces
    E_dsdna: float = 250  # Young's Modulus of dsDNA is ~ 250 pn / nm^2 (0.25 GPa )
    C_repulsion: float = 2.0

    # Pre-Initialized Class Variables that users should NOT change as they are manipulated during Shape Annealing:
    node_count: int = 0
    dna_pitch_per_rise: float = 0.342  # BDNA is .342 nm in rise per basepair of length
    smallest_edge: tuple = None  # initialize a variable if we are using these for tracking later.
    NPs: list = field(default_factory=list)
    all_faces: dict = field(default_factory=dict)
    plotting_faces: dict = field(default_factory=dict)  # Will be used for plotting and distinguishing faces
    debug_mode: bool = False  # Used during debugging to check to see what issues might be happening
    terminal_edges: list = field(default_factory=list)  # List of all terminal edges so we know not to change these.

    # POST-INITIALIZATION METHOD THAT CREATES OUR INITIAL GRAPH!!!
    def __post_init__(self):
        # Verify the routing software validity:
        if self.routing_software not in ['DAEDALUS', 'vHelix']:
            raise Exception(f'INVALID Routing Software ({self.routing_software}) choice, we only support DAEDALUS, '
                            f'vHelix!')
        if self.routing_software == 'DAEDALUS':
            if self.min_edge_length not in [31, 38, 42, 52, 63, 73, 84, 94, 105, 115, 126]:
                raise Exception(f'Invalid minimal edge length. Since you are using DAEDALUS, you need to have'
                                f'a minimal edge length in units of basepairs supported by these softwares. The valid'
                                f'values are: {[42, 52, 63, 73, 84, 94, 105, 115, 126]}')
        elif self.routing_software == 'vHelix':
            if self.min_edge_length < 8:
                raise Exception(f'Invalid minimal edge length of {self.min_edge_length}. The minimal edge length for '
                                f'vHelix is 8nm.')

        self.scaling_edge = None  # This is an edge that is used if DAEDALUS is used.
        self.design_graph = nx.empty_graph()  # Initialize this to an empty graph
        # Initialize the plotting-faces dictionary:
        self.plotting_faces['NP'], self.plotting_faces['binding_region'], self.plotting_faces['other'] = [], [], []
        self.valid_to_merge = []  # List used to check which faces are eligible to merge
        self.nonterminal_nodes = []
        self.edge_divide_nodes = {}
        self.max_edge_length = self.min_edge_length * self.max_edge_length_multiplier
        self.comparison_graph = nx.empty_graph()  # Initialize an empty graph for use in objective function

        # Define out the "bounding region" as a separate graph to plot:
        self.bounding_graph = nx.empty_graph()  # Initialize as empty
        nodes = [('b1', {'bounding': True, 'x': 0, 'y': 0, 'z': 0}),
                 ('b2', {'bounding': True, 'x': self.maxX, 'y': 0, 'z': 0}),
                 ('b3', {'bounding': True, 'x': self.maxX, 'y': self.maxY, 'z': 0}),
                 ('b4', {'bounding': True, 'x': 0, 'y': self.maxY, 'z': 0}),
                 ('b5', {'bounding': True, 'x': 0, 'y': 0, 'z': self.maxZ}),
                 ('b6', {'bounding': True, 'x': self.maxX, 'y': 0, 'z': self.maxZ}),
                 ('b7', {'bounding': True, 'x': self.maxX, 'y': self.maxY, 'z': self.maxZ}),
                 ('b8', {'bounding': True, 'x': 0, 'y': self.maxY, 'z': self.maxZ})
                 ]
        self.bounding_graph.add_nodes_from(nodes)
        edges = [('b1', 'b2'), ('b2', 'b3'), ('b3', 'b4'), ('b4', 'b1'), ('b1', 'b5'), ('b4', 'b8'), ('b3', 'b7'),
                 ('b2', 'b6'), ('b5', 'b6'), ('b6', 'b7'), ('b7', 'b8'), ('b8', 'b5')]
        self.bounding_graph.add_edges_from(edges)

        # Now we call the function to create the NP spaces:
        all_pts = self.create_NP_spaces()
        if all_pts is not None:
            self.og_NP_IDs = []
            for i in all_pts.values():
                self.og_NP_IDs.append(i)
        self.create_binding_regions(all_pts=all_pts)

        # Next, we want to make a note of all the vertices and edges used. Since this is a variable thing, doing it in
        # the post_init makes the most sense:
        self.initial_nodes_pre_start = []
        for n in self.design_graph.nodes:
            self.initial_nodes_pre_start.append(n)

        # Print our error or warnings of things that aren't working yet:
        self.print_errors_or_things_im_aware_of()

        # Next verify the geometry:
        if not self.debug_mode:  # Temporary hopefully i don't forget to remove this.
            self.verify_initial_geometry()

    ####################################################################################################################
    ####################################### Functions used in Post Init ################################################
    ####################################################################################################################
    def create_NP_spaces(self) -> dict:
        """
        This function will use the passed in NP_Regions and create the spaces for these in our 3D Graph
        :return: It will return a dict of points used in the NP creation, this will become important in the binding
                 region initialization
        """
        # First check to make sure nanoparticles were properly passed in:
        if len(self.NP_Regions) == 0:
            return
            #raise ValueError('No passed in preserved geometries so there is nothing to optimize for, CANCELLING RUN!')

        # Next we loop over all NPs in the list and construct variables within our master design graph:
        all_pts = {}
        for i in self.NP_Regions:
            numNodes = i.NP_as_graph.number_of_nodes()
            # We loop over the number of nodes in the region and change the relabel dict values as needed:
            for j in range(0, numNodes):
                pt = (i.graph_points[j]['x'], i.graph_points[j]['y'], i.graph_points[j]['z'])
                if pt not in all_pts:
                    all_pts[pt] = self.node_count
                    i.relabel_dict[j] = self.node_count
                    self.increment_node_count()
                else:
                    # If the point already exists, we simply change the relabel dict to this value:
                    i.relabel_dict[j] = all_pts[pt]


            # Now we can create a new graph which is the relabelled version of the current i value of preserved geom.
            new_NP = nx.relabel_nodes(i.NP_as_graph, i.relabel_dict)
            # Finally we just add these new nodes and edges in from new_NP to our design_graph. We add the nodes in
            # the following fashion so that their attribute values carry over (since our attributes do match)
            self.design_graph.add_nodes_from([(n, new_NP.nodes[n]) for n in new_NP.nodes])
            self.design_graph.add_edges_from(new_NP.edges)

            # Add this NP's vertices to our master list to use later:
            for key, face in i.faces.items():
                # First, we need to use the relabel_dict to rename the values in the key:
                new_key = []
                for k in list(key):
                    new_key.append(i.relabel_dict[k])
                # Now we just append to our faces list:
                self.all_faces[tuple(new_key)] = face
                self.plotting_faces['NP'].append(tuple(new_key))  # Add the faces for the NP to plot
        return all_pts

    def create_binding_regions(self, all_pts: dict) -> None:
        """
        This function will used the passed in binding region spaces and specify these regions in our plot
        :param all_pts: These are the dictionary of points used in the NP Region initializaiton
        :return: Nothing this time, just updating our master graph
        """
        if all_pts is None:
            all_pts = {}

        for i in self.binding_regions:
            numNodes = i.binder_as_graph.number_of_nodes()
            # We loop over the number of nodes in the region and change the relabel dict values as needed:
            for j in range(0, numNodes):
                pt = (i.graph_points[j]['x'], i.graph_points[j]['y'], i.graph_points[j]['z'])
                if pt not in all_pts:
                    all_pts[pt] = self.node_count
                    i.relabel_dict[j] = self.node_count
                    self.increment_node_count()
                else:
                    # If the point already exists, we simply change the relabel dict to this value:
                    i.relabel_dict[j] = all_pts[pt]
            # Merging graphs, documentation of these functions found in create_NP_spaces as it is similar code.
            new_binder = nx.relabel_nodes(i.binder_as_graph, i.relabel_dict)
            self.design_graph.add_nodes_from([(n, new_binder.nodes[n]) for n in new_binder.nodes])
            self.design_graph.add_edges_from(new_binder.edges)

            # Now, if this was a face (triangular or square) we also add to our faces list:
            if i.face is not None:
                new_key = []
                for k in list(i.face[0]):
                    new_key.append(i.relabel_dict[k])
                self.all_faces[tuple(new_key)] = i.face[1]
                self.plotting_faces['binding_region'].append(tuple(new_key))

    def initialize_random_seed(self, rng_seed: int) -> None:
        """
        Simply initialize the seed for the graph:
        :return:
        """
        #if rng_seed is None:
            #rng_seed = self.random_seed
        seed(rng_seed)
        np.random.seed(rng_seed)
        np.set_printoptions(precision=3, suppress=True)

    def verify_initial_geometry(self) -> None:
        """
        This function is responsible for checking the initial geometry to verify that we can perform optimization on it.
        For example, if there is an edge length input that is < min_edge_length THEN we have an issue!
        :return:
        """
        c1 = self.outside_design_space()
        c2 = self.vert_in_NP_zone()
        c3 = self.edge_in_NP_zone()
        c4 = self.invalid_edge_lengths()
        c5 = self.used_too_much_scaffold()
        c6 = self.invalid_face_angles()
        c7 = self.too_many_edges_at_a_node()

        if c1:
            if self.debug_mode:
                self.display_graph()
            raise Exception('ERROR: You have specified geometries outside the max X, Y, or Z values.')
        elif c2:
            if self.debug_mode:
                self.display_graph()
            raise Exception('ERROR: You have specified geometries whose vertices are inside of Nanoparticle Regions.')
        elif c3:
            if self.debug_mode:
                self.display_graph()
            raise Exception('ERROR: You have specified geometries whose edges are inside of Nanoparticle Regions.')
        elif c4:
            if self.debug_mode:
                self.display_graph()
            raise Exception('ERROR: You have specified geometries whose edges are too short or too long. If the '
                            'geometry you specified does not contain small edges, it is probably due to the Convex '
                            'Hull start shape causing an issue; try simplifying your design or expanding the design '
                            'space.')
        elif c5:
            if self.debug_mode:
                self.display_graph()
            raise Exception('ERROR: You have specified geometries whose edges are too short or too long!')

        elif c6:
            if self.debug_mode:
                self.display_graph()
            raise Exception(f'ERROR: The Cnovex Hull start shape has invalid face angles. Either reduce the value of '
                            f'min_face_angle (your is set to {self.min_face_angle}), or change the starting geometry!')

        elif c7:
            if self.debug_mode:
                self.display_graph()
            raise Exception(f'ERROR: The start shape is not able to create a design with a number of edges less than '
                            f'the maximum value of: {self.max_number_edges_at_node}. You can increase this value, but '
                            f'if it is greater than 6 we have seen reduced effectiveness in generated designs. ')

    @staticmethod
    def print_errors_or_things_im_aware_of() -> None:
        """
        This function should be deleted (DeleteMe, delete me, Delete Me, delete me) once all of its components are
        fixed up!!!
        :return:
        """
        print('ERRORS OR ISSUES I AM AWARE OF:')
        print('Add a logging module AJ!')

    ####################################################################################################################
    ########################################### General Class Methods ##################################################
    ####################################################################################################################
    def increment_node_count(self):
        """
        Increments the node_count value by 1 so we can track the node values
        :return: None
        """
        self.node_count += 1


    def add_new_node(self, x: float, y: float, z: float, terminal: bool) -> None:
        """
        This function will create a new node in our design graph
        :param x: X position
        :param y: Y Position
        :param z: Z position
        :param terminal: Is node terminal or not
        :return: Nothing, just adds a new node.
        """
        self.design_graph.add_node(self.node_count, x=x, y=y, z=z, terminal=terminal)
        self.increment_node_count()


    def display_graph(self, display_node_name: bool = False, only_show_start_criteria: bool = False,
                      save_or_show: str = 'show', savepath: str = None, return_data: bool = False,
                      hide_background: bool = False, return_figure: bool = False, use_in_dash: bool = False):
        """
        This method will take the current nodes and edges and plot them so we can visualize the current design space
        :param display_node_name: Display the node numbers as a text value based on their value.
        :return: None
        """
        if save_or_show == 'save' and savepath is None:
            warnings.warn('No savepath was specified, the graph will be shown instead and you can manually save it.')
        fig = go.Figure()
        # First displaying the bounding cube region based on user defined inputs:
        first_plot = True

        # Determining if we are showing text or not on the points:
        if display_node_name:
            mode = "markers+text"
        else:
            mode = "markers"
        # First plotting the nodes of the bounding cube:
        for node in self.bounding_graph.nodes():
            n = self.bounding_graph.nodes[node]
            x, y, z = n['x'], n['y'], n['z']
            # Creating the traces:
            if first_plot:
                node_trace = go.Scatter3d(x=[x], y=[y], z=[z], mode=mode, text=f'{node}',
                                          marker=dict(size=6, color='black', symbol='square'),
                                          textposition='top center', name='Bounding Region',
                                          legendgroup='Bounding Region', textfont=dict(size=12))
                first_plot = False
            else:
                node_trace = go.Scatter3d(x=[x], y=[y], z=[z], mode=mode, text=f'{node}',
                                          marker=dict(size=6, color='black', symbol='square'),
                                          textfont=dict(size=12), name='Bounding Region', textposition='top center',
                                          legendgroup='Bounding Region', showlegend=False)
            # Add to the plot:
            fig.add_trace(node_trace)

        # Next the edges of the bounding cube:
        for edge in self.bounding_graph.edges():
            n1, n2 = self.bounding_graph.nodes[edge[0]], self.bounding_graph.nodes[edge[1]]
            x1, y1, z1 = n1['x'], n1['y'], n1['z']
            x2, y2, z2 = n2['x'], n2['y'], n2['z']
            node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                      line=dict(color='black', dash='dash', width=6),
                                      name='Bounding Region', legendgroup='Bounding Region', showlegend=False)
            fig.add_trace(node_trace)

        # Next up we want to plot the edges found within the Design Graph:
        first_edge_plot = True
        first_binding_edge_plot = True
        first_mixed_edge = True
        first_nonterminal_edge = True

        # For loop used in the publication print of these figures to find binding edges and things:
        edge_node_pairs = []
        initial_plot_edges = []  # Used to visualize the NP space when showing the pre-triangulation criteria
        for i in self.binding_regions:
            if i.binding_type == 'edge':
                p1, p2 = self.node_from_xyz(xyz=i.v1), self.node_from_xyz(xyz=i.v2)
                edge_node_pairs.extend([(p1, p2), (p2, p1)])
            if i.binding_type == 'triangle' or i.binding_type == 'tri':
                p1, p2, p3 = self.node_from_xyz(xyz=i.v1), self.node_from_xyz(xyz=i.v2), self.node_from_xyz(xyz=i.v3)
                edge_node_pairs.extend([(p1, p2), (p2, p1), (p1, p3), (p3, p1), (p3, p2), (p2, p3)])
            if i.binding_type == 'rect':
                p1, p2 = self.node_from_xyz(xyz=i.v1), self.node_from_xyz(xyz=i.v2)
                p3, p4 = self.node_from_xyz(xyz=i.v3), self.node_from_xyz(xyz=i.v4)
                edge_node_pairs.extend([(p1, p2), (p2, p1), (p1, p4), (p4, p1), (p3, p2), (p2, p3), (p3, p4), (p4, p3)])
        if self.plotting_faces['other']:  # Only showing green lines when a polyhedra is generated
            for i in self.plotting_faces['NP']:
                edge_node_pairs.extend([(i[0], i[1]), (i[1], i[0]), (i[0], i[3]), (i[3], i[0]), (i[2], i[1]),
                                        (i[1], i[2]), (i[2], i[3]), (i[3], i[2])])
        for i in self.plotting_faces['NP']:
            initial_plot_edges.extend([(i[0], i[1]), (i[1], i[0]), (i[0], i[3]), (i[3], i[0]), (i[2], i[1]),
                                       (i[1], i[2]), (i[2], i[3]), (i[3], i[2])])


        for edge in self.design_graph.edges():
            n1, n2 = self.design_graph.nodes[edge[0]], self.design_graph.nodes[edge[1]]
            x1, y1, z1 = n1['x'], n1['y'], n1['z']
            x2, y2, z2 = n2['x'], n2['y'], n2['z']
            if self.debug_mode:
                if first_edge_plot and n1['terminal'] and n2['terminal'] and 'binding_zone' not in n1 and \
                        'binding_zone' not in n2:
                    node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                              line=dict(color='black', width=4),
                                              name='Graph Edges (Both Terminal)', legendgroup='Graph Edges (Both Terminal)')
                    first_edge_plot = False
                elif not first_edge_plot and n1['terminal'] and n2['terminal'] and 'binding_zone' not in n1 and \
                        'binding_zone' not in n2:
                    node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                              line=dict(color='black', width=4), showlegend=False,
                                              name='Graph Edges (Both Terminal)', legendgroup='Graph Edges (Both Terminal)')
                elif first_binding_edge_plot and 'binding_zone' in n1 and 'binding_zone' in n2:
                    node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                              line=dict(color='blue', width=4),
                                              name='Graph Edges (Both Binding Nodes)', legendgroup='Binding Edges')
                    first_binding_edge_plot = False
                elif not first_binding_edge_plot and 'binding_zone' in n1 and 'binding_zone' in n2:
                    node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                              line=dict(color='blue', width=4), showlegend=False,
                                              name='Graph Edges (Both Binding Nodes)', legendgroup='Binding Edges')
                elif first_nonterminal_edge and not n1['terminal'] and not n2['terminal']:
                    node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                              line=dict(color='cyan', width=4),
                                              name='Graph Edges (Both Non-Terminal Nodes)', legendgroup='NonT Edges')
                    first_nonterminal_edge = False
                elif not first_nonterminal_edge and not n1['terminal'] and not n2['terminal']:
                    node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                              line=dict(color='cyan', width=4), showlegend=False,
                                              name='Graph Edges (Both Non-Terminal Nodes)', legendgroup='NonT Edges')
                else:
                    ## This last case is when one of the nodes is a binding point and the other is not, so just showing
                    ## this normally for now with a black edge for a mixed edge:
                    node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                              line=dict(color='yellow', width=4), showlegend=first_mixed_edge,
                                              name='Graph Edges (Mixed Nodes)', legendgroup='Mixed Edge')
                    first_mixed_edge = False
            else:
                if edge in edge_node_pairs and edge not in initial_plot_edges:
                    node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                              line=dict(color='green', width=12), showlegend=first_edge_plot,
                                              name='Graph Edge', legendgroup='Graph Edge')
                    first_edge_plot = False
                elif edge in initial_plot_edges:
                    node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                              line=dict(color='red', width=12), showlegend=first_edge_plot,
                                              name='Graph Edge', legendgroup='Graph Edge')
                    first_edge_plot = False
                else:
                    node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                              line=dict(color='black', width=12), showlegend=first_edge_plot,
                                              name='Graph Edge', legendgroup='Graph Edge')
                    first_edge_plot = False

            fig.add_trace(node_trace)

        # Plotting the vertices found within the Design Graph:
        first_non_t_node_plot = True
        first_terminal_node = True
        first_binding_node = True
        first_trace = True

        # For loop used for finding points in the NPs (if any NPs)
        verts_in_NPs = []
        for n in self.plotting_faces['NP']:
            for i in list(n):
                if i not in verts_in_NPs:
                    verts_in_NPs.append(i)
        for node in self.design_graph.nodes():
            n = self.design_graph.nodes[node]
            if self.debug_mode:
                if first_non_t_node_plot and not self.design_graph.nodes[node]['terminal'] and 'binding_zone' not in \
                        self.design_graph.nodes[node]:
                    node_trace = go.Scatter3d(x=[n['x']], y=[n['y']], z=[n['z']], mode=mode, text=f'{node}',
                                              marker=dict(size=6, color='fuchsia'), name="Non-Terminal Nodes",
                                              legendgroup="Non-Terminal Nodes")
                    first_non_t_node_plot = False
                elif not first_non_t_node_plot and not self.design_graph.nodes[node]['terminal'] and 'binding_zone' \
                        not in self.design_graph.nodes[node]:
                    node_trace = go.Scatter3d(x=[n['x']], y=[n['y']], z=[n['z']], mode=mode, text=f'{node}',
                                              marker=dict(size=6, color='fuchsia'), name="Non-Terminal Nodes",
                                              legendgroup="Non-Terminal Nodes", showlegend=False)

                elif first_terminal_node and 'terminal' in self.design_graph.nodes[node] and 'binding_zone' not in \
                        self.design_graph.nodes[node]:
                    node_trace = go.Scatter3d(x=[n['x']], y=[n['y']], z=[n['z']], mode=mode, text=f'{node}',
                                              marker=dict(size=6, color='red'), name="Terminal Nodes",
                                              legendgroup="Terminal Nodes")
                    first_terminal_node = False
                elif not first_terminal_node and 'terminal' in self.design_graph.nodes[node] and 'binding_zone' not in \
                        self.design_graph.nodes[node]:
                    node_trace = go.Scatter3d(x=[n['x']], y=[n['y']], z=[n['z']], mode=mode, text=f'{node}',
                                              marker=dict(size=6, color='red'), name="Terminal Nodes",
                                              legendgroup="Terminal Nodes", showlegend=False)
                elif first_binding_node and 'binding_zone' in self.design_graph.nodes[node]:
                    node_trace = go.Scatter3d(x=[n['x']], y=[n['y']], z=[n['z']], mode=mode, text=f'{node}',
                                              marker=dict(size=6, color='orange', symbol='square'), name="Binding Point(s)",
                                              legendgroup="Binding Point(s)")
                    first_binding_node = False
                else:
                    node_trace = go.Scatter3d(x=[n['x']], y=[n['y']], z=[n['z']], mode=mode, text=f'{node}',
                                              marker=dict(size=6, color='orange', symbol='square'), name="Binding Point(s)",
                                              legendgroup="Binding Point(s)", showlegend=False)
            else:
                if node in verts_in_NPs:
                    node_trace = go.Scatter3d(x=[n['x']], y=[n['y']], z=[n['z']], mode=mode, text=f'{node}',
                                              marker=dict(size=12, color='red', symbol='square'),
                                              name="Terminal Nodes",
                                              legendgroup="Terminal Nodes", showlegend=first_trace)
                    first_trace = False
                else:
                    node_trace = go.Scatter3d(x=[n['x']], y=[n['y']], z=[n['z']], mode=mode, text=f'{node}',
                                              marker=dict(size=12, color='green', symbol='square'), name="Terminal Nodes",
                                              legendgroup="Terminal Nodes", showlegend=first_trace)
                    first_trace = False


            fig.add_trace(node_trace)


        # Plotting the Nanoparticle Faces:
        for nap in self.NP_Regions:
            for f in nap.faces.values():
                nt = f.get_face_plot(color_type='NP', override_opacity=True)
                fig.add_trace(nt)


        # Plotting the faces:
        if not only_show_start_criteria:
            for face_type, face_verts in self.plotting_faces.items():
                if face_type == 'NP':
                    pass  # I already plot these above.
                else:
                    # We need to loop through and break down each face individual to plot it:
                    for face in face_verts:
                        if len(face) == 4:
                            plot_face = MeshFaces(v1=self.xyz_from_node(node=face[0]), v2=self.xyz_from_node(node=face[1]),
                                                  v3=self.xyz_from_node(node=face[2]), v4=self.xyz_from_node(node=face[3]),
                                                  face_type='rect', faces=face)
                        else:
                            plot_face = MeshFaces(v1=self.xyz_from_node(node=face[0]), v2=self.xyz_from_node(node=face[1]),
                                                  v3=self.xyz_from_node(node=face[2]), face_type='triangle')

                        # Now, we simply add a new trace based on the mesh face:
                        new_trace = plot_face.get_face_plot(color_type=face_type)
                        fig.add_trace(new_trace)


            # Since we use Mesh3D plots, we need to create "dummy" plts for the legend to show what is what:
            dummy_scatter1 = go.Scatter3d(x=[None], y=[None], z=[None], name='Nanoparticle Faces', mode='markers',
                                          marker=dict(color='red', opacity=0.4, size=8, symbol='square'))
            dummy_scatter2 = go.Scatter3d(x=[None], y=[None], z=[None], name='Binding Region Faces', mode='markers',
                                          marker=dict(color='green', opacity=0.4, size=8, symbol='square'))
            dummy_scatter3 = go.Scatter3d(x=[None], y=[None], z=[None], name='Output Polygon', mode='markers',
                                          marker=dict(color='grey', opacity=0.4, size=8, symbol='square'))
            fig.add_trace(dummy_scatter1)
            fig.add_trace(dummy_scatter2)
            fig.add_trace(dummy_scatter3)
        # Lastly, update the axis titles and their values:
        global axis_font_size
        global title_font_size

        camera = dict(
            eye=dict(x=-2, y=-0.26, z=0.6)
        )

        # Update the layout with axis titles and annotations
        if display_node_name:
            fig.update_layout(
                title=dict(text="Graph Representation of PLY File Generated via Shape Annealing",
                           font=dict(size=title_font_size)),
                font=dict(size=axis_font_size),
                scene=dict(xaxis=dict(title="X Axis (nm)"),
                           yaxis=dict(title="Y Axis (nm)"),
                           zaxis=dict(title="Z Axis (nm)")),
                # annotations=annotations,
                scene_camera=camera
            )
        else:
            fig.update_layout(
                title=dict(text="Graph Representation of PLY File Generated via Shape Annealing",
                           font=dict(size=title_font_size)),
                scene=dict(
                    xaxis=dict(title="X (nm)", title_font=dict(size=axis_font_size+4), tickfont=dict(size=18)),
                    yaxis=dict(title="Y (nm)", title_font=dict(size=axis_font_size+4), tickfont=dict(size=18)),
                    zaxis=dict(title="Z (nm)", title_font=dict(size=axis_font_size+4), tickfont=dict(size=18))
                ),
                scene_camera=camera
            )

        if hide_background:
            fig.update_layout(
                scene=dict(
                    xaxis_visible=False,
                    yaxis_visible=False,
                    zaxis_visible=False,
                    bgcolor="white",  # Set background color to white
                ),
                plot_bgcolor="white",  # Set plot background color to white
                paper_bgcolor="white"  # Set paper background color to white
            )

        if use_in_dash:
            camera = dict(
                eye=dict(x=-2, y=-0.26, z=0.6)
            )
            fig.update_layout(
                title=None,
                font=dict(size=18),
                scene=dict(xaxis=dict(title="X Axis (nm)", showbackground=False),
                           yaxis=dict(title="Y Axis (nm)", showbackground=False),
                           zaxis=dict(title="Z Axis (nm)", showbackground=False),
                           bgcolor="white"),
                plot_bgcolor="white",  # Set plot background color to white
                paper_bgcolor="white",  # Set paper background color to white
                showlegend=False,
                scene_camera=camera
            )

        # If we are returning data for the design evolution, we just return the figure:
        if return_data:
            data_to_return = fig.data
            return data_to_return

        if return_figure:
            return fig

        # Show the figure:
        if save_or_show == 'show':
            fig.show()
        elif save_or_show == 'save' and savepath is None:
            fig.show()
        else:
            # Save the interactive HTML File:
            HTML_File = os.path.join(savepath, '3D_Graph_Plot.html')
            pyo.plot(fig, filename=HTML_File, auto_open=False)
            # Save a still image:
            image_file = os.path.join(savepath, '3D_Graph_Plot.png')
            pio.write_image(fig, image_file, format="png")

    def display_graph_as_cylindrical_rep(self, hide_background: bool = True, return_figure: bool = False,
                                         hide_title: bool = True):
        CYLINDER_RESOLUTION = 25
        def rotation_matrix_from_vectors(vec1, vec2):
            a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
            v = np.cross(a, b)
            s = np.linalg.norm(v)
            if s < 1e-8:
                # When vectors are almost parallel, align the cylinder with the direction vector
                rotation_matrix = np.eye(3)
            else:
                # Calculate the rotation matrix
                kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
                rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - np.dot(a, b)) / (s ** 2))
            return rotation_matrix

        fig = go.Figure()
        # First plotting the nodes of the bounding cube:
        mode = "markers"
        for node in self.bounding_graph.nodes():
            n = self.bounding_graph.nodes[node]
            x, y, z = n['x'], n['y'], n['z']
            # Creating the traces:
            node_trace = go.Scatter3d(x=[x], y=[y], z=[z], mode=mode, text=f'{node}',
                                      marker=dict(size=6, color='black', symbol='square'),
                                      textfont=dict(size=12), name='Bounding Region', textposition='top center',
                                      legendgroup='Bounding Region', showlegend=False)
            # Add to the plot:
            fig.add_trace(node_trace)

        # Next the edges of the bounding cube:
        for edge in self.bounding_graph.edges():
            n1, n2 = self.bounding_graph.nodes[edge[0]], self.bounding_graph.nodes[edge[1]]
            x1, y1, z1 = n1['x'], n1['y'], n1['z']
            x2, y2, z2 = n2['x'], n2['y'], n2['z']
            node_trace = go.Scatter3d(x=[x1, x2], y=[y1, y2], z=[z1, z2], mode='lines',
                                      line=dict(color='black', dash='dash', width=6),
                                      name='Bounding Region', legendgroup='Bounding Region', showlegend=False)
            fig.add_trace(node_trace)

        # Next up we want to plot the edges found within the Design Graph:
        r = self.dna_diameter / 2
        for edge in self.design_graph.edges():
            n1, n2 = self.design_graph.nodes[edge[0]], self.design_graph.nodes[edge[1]]
            P1 = np.array([n1['x'], n1['y'], n1['z']])
            P2 = np.array([n2['x'], n2['y'], n2['z']])
            direction = P2 - P1
            length = np.linalg.norm(direction)

            rotation_matrix = rotation_matrix_from_vectors(np.array([0, 0, 1]), direction)
            theta = np.linspace(0, 2 * np.pi, CYLINDER_RESOLUTION)
            z = np.linspace(0, length, CYLINDER_RESOLUTION)
            theta_grid, z_grid = np.meshgrid(theta, z)
            cylinder_surface_points = np.array([r * np.cos(theta_grid), r * np.sin(theta_grid), z_grid])
            if np.allclose(rotation_matrix, np.eye(3)):
                # Handle the case when vectors are almost parallel
                pass
            else:
                # Apply the rotation matrix to cylinder surface points
                cylinder_surface_points = np.array([r * np.cos(theta_grid), r * np.sin(theta_grid), z_grid])
                for i in range(CYLINDER_RESOLUTION):
                    for j in range(CYLINDER_RESOLUTION):
                        cylinder_surface_points[:, i, j] = np.dot(rotation_matrix, cylinder_surface_points[:, i, j])

            # Translate the points to the desired position in 3D space
            cylinder_surface_points = cylinder_surface_points + P1.reshape(-1, 1, 1)
            x_grid, y_grid, z_grid = cylinder_surface_points[0], cylinder_surface_points[1], cylinder_surface_points[2]
            x_points = x_grid.reshape(-1)
            y_points = y_grid.reshape(-1)
            z_points = z_grid.reshape(-1)

            # Create vertices for the cylinder's surface
            faces, vertices = [], []
            for i in range(CYLINDER_RESOLUTION):
                for j in range(CYLINDER_RESOLUTION):
                    x = x_points[i * CYLINDER_RESOLUTION + j]
                    y = y_points[i * CYLINDER_RESOLUTION + j]
                    z = z_points[i * CYLINDER_RESOLUTION + j]
                    vertices.append([x, y, z])

            # Create faces by defining the vertex indices
            for i in range(CYLINDER_RESOLUTION-1):
                for j in range(CYLINDER_RESOLUTION-1):
                    v1 = i * CYLINDER_RESOLUTION + j
                    v2 = v1 + 1
                    v3 = (i + 1) * CYLINDER_RESOLUTION + j
                    v4 = v3 + 1
                    faces.extend([(v1, v2, v3), (v2, v4, v3)])

            fig.add_trace(go.Mesh3d(x=x_points,
                                    y=y_points,
                                    z=z_points,
                                    i=[face[0] for face in faces], j=[face[1] for face in faces],
                                    k=[face[2] for face in faces], color='darkgray', showlegend=False))


        # Plotting the Nanoparticle Faces:
        for nap in self.NP_Regions:
            for f in nap.faces.values():
                nt = f.get_face_plot(color_type='NP', override_opacity=True)
                fig.add_trace(nt)

        if hide_background:
            if hide_title:
                title_text = None
            else:
                title_text = dict(text="Cylindrical Representation of DNA Origami",
                                  font=dict(size=title_font_size))
            camera = dict(
                eye=dict(x=-2, y=-0.26, z=0.6)
            )
            fig.update_layout(
                title=title_text,
                font=dict(size=18),
                scene=dict(xaxis=dict(title="X Axis (nm)", showbackground=False),
                           yaxis=dict(title="Y Axis (nm)", showbackground=False),
                           zaxis=dict(title="Z Axis (nm)", showbackground=False),
                           bgcolor="white"),
                plot_bgcolor="white",  # Set plot background color to white
                paper_bgcolor="white",  # Set paper background color to white
                scene_camera=camera
            )
        if return_figure:
            return fig
        else:
            fig.show()

    def get_all_verts(self, filter: list = None) -> list:
        """
        This function returns all the vertices of the graph in X Y Z Format
        :param filter: This will be a list of verts to retrive
        :return: list of points
        """
        all_verts = []
        if filter is not None:
            # If we have a filter of points, we only get the points in the passed in list.
            unique_digits = set()
            for p in filter:
                unique_digits.update(p)
            # Loop over the unique digits and retrieve these points:
            for digit in list(unique_digits):
                pt = self.xyz_from_node(node=digit)
                all_verts.append(pt)
        else:
            # Otherwise we just grab all the verts to return
            for node in self.design_graph.nodes():
                pt = self.xyz_from_node(node=node)
                all_verts.append(pt)
        return all_verts

    def xyz_from_node(self, node: int) -> np.array:
        """
        This function will take in a node number and return that points X, Y, Z coordinates as a numpy array
        :param node: Node number we want the X Y Z Values of
        :return: [X Y Z]
        """
        return np.array([self.design_graph.nodes[node]['x'], self.design_graph.nodes[node]['y'],
                         self.design_graph.nodes[node]['z']])

    def node_from_xyz(self, xyz: np.array) -> int:
        """
        This function takes in the coordinates of X Y Z and returns the node number from the design graph
        :param xyz: NumPy array in form of [X, Y, Z]
        :return: Integer of node number
        """
        for n in self.design_graph.nodes:
            if self.design_graph.nodes[n]['x'] == xyz[0] and self.design_graph.nodes[n]['y'] == xyz[1] and \
                    self.design_graph.nodes[n]['z'] == xyz[2]:
                return n
        raise Exception('Unable to find node in the design graph!!!')
        pass

    def nodes_on_same_face(self, n1: int, n2: int) -> bool:
        """
        This function will take in two nodes and return TRUE if they ARE on the same face, and FALSE if not
        :param n1: 1st node to check
        :param n2: 2nd node to check
        :return: True: Nodes on same face; False: Nodes NOT on same face.
        """
        for f in self.all_faces.keys():
            if n1 in list(f) and n2 in list(f):
                # If both nodes are in a tuple that is defining a face, then the edge must exist.
                return True
        # If both nodes are not in any face, then we have a new edge.
        return False

    def calculate_edge_length(self, n1: int, n2: int) -> float:
        """
        This function strictly calculated the distance between two passed in nodal points from our graph.
        :param n1: Node 1 point
        :param n2: Node 2 point
        :return: Length of the edge
        """
        global numDecimals
        p1 = [self.design_graph.nodes[n1]['x'], self.design_graph.nodes[n1]['y'], self.design_graph.nodes[n1]['z']]
        p2 = [self.design_graph.nodes[n2]['x'], self.design_graph.nodes[n2]['y'], self.design_graph.nodes[n2]['z']]
        return np.round(np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2), numDecimals)

    def all_faces_from_node(self, node: int) -> list:
        """
        This function will take in a node number and return a list of all the faces that use this node
        :param node:
        :return:
        """
        faces_with_node = []
        all_faces = self.plotting_faces['other']  # We can not merge the NP or Binding Region faces so just these
        for i in all_faces:
            if node in list(i):
                faces_with_node.append(i)
        return faces_with_node

    def group_coplanar_faces(self, faces: list, ref_pt: np.array):
        """
        This function will take a list of faces found at a specific point and return lists of any coplanar faces
        :param ref_pt: Reference point to use in planar calcs
        :param faces: Faces found at some node
        :return: IDK yet
        """
        face_groups = defaultdict(list)
        for face in faces:
            f = self.all_faces[face]
            face_norm = f.review_vertex_orientation(reference_point=ref_pt)
            face_groups[tuple(face_norm)].append(face)

        coplanar_faces = []
        for group in face_groups.values():
            if len(group) > 1:
                coplanar_faces.extend(group)
        return coplanar_faces

    def find_coincident_faces(self, face: tuple) -> list:
        """
        This function will take in a face as a tuple and search for the faces that are touching this face.
        :param face: Face we are looking for
        :return: list of faces that are coincident (share an edge)
        """
        coincident_faces = []
        all_faces_to_search = self.plotting_faces['other']
        for v1, v2 in itertools.combinations(list(face), 2):
            # We search for all combinations of 2 since we have triangular faces and look for faces that are touching:
            for all_other_faces in all_faces_to_search:
                if v1 in list(all_other_faces) and v2 in list(all_other_faces) and all_other_faces != face:
                    coincident_faces.append(all_other_faces)

        # Returning the coincident faces:
        return coincident_faces

    def find_angle_between_faces(self, face1: tuple, face2: tuple) -> np.array:
        """
        This function will find the angle between two coincident faces for the purpose of finding coplanar planes.
        :param face1: Tuple of face 1 vertices
        :param face2: Tuple of face 2 vertices
        :return: Angle (in degrees or radians, whatever) between the two normals of the planes
        """
        f1, f2 = self.all_faces[face1], self.all_faces[face2]
        f1_abcd, f2_abcd = f1.equation_of_face(), f2.equation_of_face()
        # I think we just need a, b and c:
        n1 = np.array([f1_abcd[0], f1_abcd[1], f1_abcd[2]])
        n2 = np.array([f2_abcd[0], f2_abcd[1], f2_abcd[2]])

        # Now, due to rounding errors in calculations (for example, if we evaluate the below to -1.0000000002, we want
        # to round the value:
        global numDecimals
        angle_of_planes = np.round(np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2)), numDecimals)
        return np.degrees(np.arccos(angle_of_planes))

    def are_faces_coincident(self, face1: list, face2: list) -> bool:
        """
        This function is in support of find_coplanar_faces to double check that when we are combining similar
        coplanar faces that the faces are also coincident in some capacity. This is because two faces may be 180
        degrees apart BUT the faces may NOT be coincident.
        :param face1: List of faces to check
        :param face2:
        :return:
        """
        unique_1 = list(set(value for tpl in face1 for value in tpl))
        unique_2 = list(set(value for tpl in face2 for value in tpl))
        common_values = set(unique_1).intersection(unique_2)
        if len(common_values) >= 2:
            # Basically, if the faces we are searching share at least 2 nodes that means they share an edge and so these
            # faces must be coincident
            return True
        else:
            # Otherwise if they only share 1 or 0, they are not coincident and just happen to be 180 degree apart planes
            # which are NOT coincident.
            return False

    def find_coplanar_faces(self) -> None:
        # Reset these values:
        self.valid_to_merge = []
        possible_merges = []
        # To be able to find coplanar faces, we need to understand face orientation, so we use a reference point which
        # i will take as the mean of all vertex values:
        all_faces_to_search = deepcopy(self.plotting_faces['other'])
        stop_while_count = 0
        # We start by running a while loop over the faces to search. We will remove faces from this last as we go:
        while all_faces_to_search and stop_while_count < 100:
            # First we select a face at random and find the faces that are coincident to any of its edges:
            curSearch = choice(all_faces_to_search)
            touching_faces = self.find_coincident_faces(face=curSearch)
            coplanar_faces = []
            # After finding the coincident faces, we want to find the angle between the two faces to see if it is
            # ~180 degrees (i.e. planar) or not
            for f in touching_faces:
                angle_between_faces = self.find_angle_between_faces(face1=curSearch, face2=f)
                if 179 <= angle_between_faces <= 181:
                    # If we are approximately 180 degrees apart, then we MUST be coplanar. I gave an error of 1 degrees.
                    coplanar_faces.extend([curSearch, f])

            # Then, after looping through here to find coplanar faces, we add to the list of possible merges:
            if coplanar_faces:
                if not possible_merges:
                    # If possible_merges is an empty list we just append
                    possible_merges.append(coplanar_faces)
                else:
                    # Otherwise, we first test if the faces we found are parallel to any others we found:
                    list_copy = deepcopy(possible_merges)
                    new_merge_list = True
                    for merge_list in list_copy:
                        angle_between_faces = self.find_angle_between_faces(face1=coplanar_faces[0],
                                                                            face2=merge_list[0])
                        # We check if the angle is about 180 degrees (with 1 degree error) meaning coplanar
                        # We also check if it is 0 since we may be comparing the normals of the same face:
                            ### NOTE: If 0 is breaking thigns AJ just make sure we don't choose two similar faces in the
                                      # above for the angle!
                        faces_coincident = self.are_faces_coincident(face1=coplanar_faces, face2=merge_list)
                        if (179 <= angle_between_faces <= 181 or np.isclose(np.round(angle_between_faces, 0), 0.)) and \
                                faces_coincident:
                            # First we note that a merge list was found so we turn off the new_merge_list:
                            new_merge_list = False
                            # If these two planes are about equivalent, then we append (but only unique values):
                            list_to_add_to = deepcopy(possible_merges[possible_merges.index(merge_list)])
                            new_list = list(set(list_to_add_to + coplanar_faces))
                            # Now reset the value for possible merges at the index:
                            possible_merges[possible_merges.index(merge_list)] = new_list
                            # And since we found a face to add to, we can break out this loop
                            break

                    if new_merge_list:
                        # And if no similar merge face was already documented, we just append the unique values because
                        # somtimes duplicates can sneak in when searching the way i am:
                        new_list = list(set(coplanar_faces))
                        possible_merges.append(new_list)

            # Then, after we search for this face and everything it is coplanar to, we just remove it remove search list
            all_faces_to_search.remove(curSearch)

            # Increment counter to prevent infinite while loop:
            stop_while_count += 1

        # After the loop, we then set the value of valid_to_merge to be our possible merge list. I get a conservative
        # answer which may have duplicates so we need to ensure they are removed:
        self.valid_to_merge = possible_merges

    def find_nonterminal_nodes(self) -> None:
        """
        This function will find all of the non-terminal nodes in the design and append to a list that will be used to
        apply a variety of production rules to manipulate geometry.
        Terminal nodes are those which belong to terminal geometries (binding regions / NP spaces). Nonterminal are
        all others
        :return: Nothing, just appends to a list
        """
        self.nonterminal_nodes = []  # Reset this list since nodes can be removed
        for n in self.design_graph.nodes:
            # If a node is not marked as terminal and not already accounted for we go into this check
            if not self.design_graph.nodes[n]['terminal'] and n not in self.nonterminal_nodes:
                self.nonterminal_nodes.append(n)

    def colinearity_check(self, v1: np.array, v2: np.array, v3: np.array) -> bool:
        """
        This function will deterine if 3 passed in points ARE colinear or not. This is useful for merging faces so we
        know which point to reconnect
        :param v1: XYZ
        :param v2: XYZ
        :param v3: XYZ
        :return: True: points ARE colinear; False: points are NOT colinear
        """
        ev1 = v2 - v1
        ev2 = v3 - v1
        cross_vector_check = np.isclose(np.cross(ev1, ev2), [0.0, 0.0, 0.0])
        if np.all(cross_vector_check):
            return True
        else:
            return False

    def update_nodal_position(self, node: int, new_xyz: np.array) -> None:
        """
        This function will take a node and it's new XYZ positions and update the NetworkX graph to represent this.
        NOTE: This does NOT update the values in the triangular faces, another function controls this.
        :param node: Node index #
        :param new_xyz: New XYZ Positions [0:X, 1:Y, 2:Z]
        :return: Nothing, just updates the design graph:
        """
        global numDecimals
        self.design_graph.nodes[node]['x'] = round(new_xyz[0], numDecimals)
        self.design_graph.nodes[node]['y'] = round(new_xyz[1], numDecimals)
        self.design_graph.nodes[node]['z'] = round(new_xyz[2], numDecimals)

    def update_face_values(self, node: int, old_xyz: np.array, new_xyz: np.array) -> None:
        """
        This function will update all of the face values in the current mesh with the new nodal values
        :param node: Node number to update
        :param old_xyz: Old XYZ Position
        :param new_xyz: New XYZ Position
        :return: Nothing, just updates values
        """
        # We loop over all the faces which contain vertex information and update these values:
        for k, v in self.all_faces.items():
            # If the node is in the list of vertices for the face, we need to update it's mesh face:
            if node in list(k):
                v.update_vertex_values(old_xyz=old_xyz, new_xyz=new_xyz)

    def remove_faces(self, face_verts: tuple) -> None:
        """
        This function removes a face from the self.all_faces and self.plotting_faces['other'].
        NOTE: this will NOT search for the NP or Binding Region faces!!!
        :param face_verts: Tuple of vertices to search for
        :return: Nothing, just updates values
        """
        self.all_faces.pop(face_verts, None)
        all_faces = self.plotting_faces['other']
        all_faces.remove(face_verts)
        self.plotting_faces['other'] = all_faces

    def get_3D_reference(self) -> np.array:
        """
        This function will return the reference point I am using to define if planes are CCW or CW orientation. This
        may need to change later on, but for now, the average point works
        :return:
        """
        vertices = np.array(self.get_all_verts())
        reference_point = np.mean(vertices, axis=0)
        return reference_point

    def calculate_3D_distance(self, p1: np.array, p2: np.array) -> float:
        """
        Simple function that is the distance formula
        :param p1: Point 1 XYZ
        :param p2: Point 2 XYZ
        :return: Distance between P1 and P2
        """
        return np.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2)

    def find_nodes_on_line(self, node: int) -> list:
        """
        This function will take in a node that was previously used as an edge divide and find all other nodes that
        belong on the same line / edge.
        :param node: Node used in an edge divide
        :return: List of nodes that are on the same edge
        """
        # First we find the other two points:
        ends = self.edge_divide_nodes[node]
        try:
            p1, p2 = self.xyz_from_node(node=ends[0]), self.xyz_from_node(node=ends[1])
        except:
            self.display_graph(display_node_name=True)
        D_AB = self.calculate_3D_distance(p1=p1, p2=p2)  # We are checking simply if these distances are equal
        points_on_line = []
        # Now we loop over the other points and find if any are on this line:
        for pt in self.design_graph.nodes:
            if pt != node:
                check_pt = self.xyz_from_node(node=pt)
                D1 = self.calculate_3D_distance(p1=p1, p2=check_pt)
                D2 = self.calculate_3D_distance(p1=p2, p2=check_pt)
                if np.isclose(D1+D2, D_AB):
                    points_on_line.append(pt)

        # At the end, we return points_on_line:
        return points_on_line

    def find_triangles_via_division(self, verts: list) -> list:
        """
        This function is going to specifically look for triangles that are touching another triangle via edge divisions
        :param verts: Vertices that belong to the same edge we are looking at
        :return: A list of faces that are touching this edge
        """
        # We just run through all combinations of 2 of these vertices. This would be computationally expensive, but
        # due to the minimal_length design constraint, there really shouldn't be that many nodes on an edge:
        faces_on_edge = []
        for v1, v2 in itertools.combinations(verts, 2):
            for i in self.all_faces.keys():
                if v1 in list(i) and v2 in list(i):
                    faces_on_edge.append(i)
        return faces_on_edge

    def calculate_triangle_area(self, n1: int, n2: int, n3: int) -> float:
        """
        This function will take the 3 vertices of a triangle and calculate it's area.
        :param n1: Node 1
        :param n2: Node 2
        :param n3: Node 3
        :return: Area of triangle
        """
        global numDecimals
        xyz1, xyz2, xyz3 = self.xyz_from_node(node=n1), self.xyz_from_node(node=n2), self.xyz_from_node(node=n3)
        AB = xyz2 - xyz1
        AC = xyz3 - xyz1
        # The magnitude of the cross product / 2 is the area of the triangle:
        mag = np.linalg.norm(np.cross(AB, AC))
        return round(float(mag / 2), numDecimals)

    def project_points(self, points: np.array) -> np.array:
        """
        THis function will utilize PCA to project 3D points into 2D using orthogonality
        :param points: List of XYZ points
        :return:
        """
        pca = PCA(n_components=2)
        return pca.fit_transform(points)

    def merge_n_faces(self, faces: list) -> None:
        """
        This function will work to merge N faces found in a plane. The logic is slightly different than merge 2 faces
        which is why this function is separated, for now at least.
        :param faces: All coplanar faces found
        :return: Nothing
        """
        # First we add all the vertices of these points to a list
        all_verts = []
        all_nodes = []
        unique_verts = list(set(value for tpl in faces for value in tpl))  # For efficiency
        for face in faces:
            for v in list(face):
                if v not in all_nodes:
                    all_verts.append(self.xyz_from_node(node=v))
                    all_nodes.append(v)
            if len(all_verts) == len(unique_verts):
                # To avoid overly long looping, we break once we find all the unique values
                break
        # Before doing a Convex Hull and Triangulation, we need to run a PCA to decompose into 2D and create a map:
        projected_verts = self.project_points(points=np.array(all_verts))
        mapping = {}
        ct = 0
        for v in projected_verts:
            mapping[tuple(np.round(v, 3))] = all_verts[ct]
            ct += 1

        # Next, we want essentially a Convex Hull that is triangulated of these points to get the MINIMAL amount of
        # triangles in return to cover some space.
        convex_hull = ConvexHull(projected_verts)
        c_hull_pts = convex_hull.points[convex_hull.vertices]
        triangulation = Delaunay(c_hull_pts)
        # And we only want the vertices of the triangles:
        triangles_to_add = triangulation.simplices

        # Now we create the arrays for our triangles:
        all_new_triangle_verts = []
        new_triangles = []
        for new_triangle in triangles_to_add:
            new_verts = []
            for new_vertex in new_triangle:
                two_d_vert = triangulation.points[new_vertex]
                # Mapping into 3D:
                three_d_vert = mapping[tuple(np.round(two_d_vert, 3))]
                new_verts.append(three_d_vert)
                # Also, to see which points need to be removed, we append to all_new_verts:
                if not any(np.all(np.isclose(three_d_vert, elem)) for elem in all_new_triangle_verts):
                    all_new_triangle_verts.append(three_d_vert)
            # After this internal loop, we append these vertices to our list tracking the new triangles:
            new_triangles.append(new_verts)

        all_plotting = self.plotting_faces['other']
        for t in new_triangles:
            node_list = []
            for v in t:
                node_list.append(self.node_from_xyz(xyz=v))
            self.all_faces[tuple(node_list)] = MeshFaces(v1=t[0], v2=t[1], v3=t[2], face_type='triangle')
            all_plotting.append(tuple(node_list))

        # Then, we need to remove the passed in faces (since they have been now merged) and remove any "lone" vertices
        # from a result of the face merging:
        for i in faces:
            all_plotting.remove(i)
            self.all_faces.pop(i, None)

        # Finding unused verts to remove:
        for vert in all_verts:
            if not any(np.array_equal(vert, new_vert) for new_vert in all_new_triangle_verts):
                node_to_remove = self.node_from_xyz(xyz=vert)
                # Before we remove this node, we first need to see if it belongs to any other triangles in other
                # faces that may not be coplanar also:
                faces_with_other_node = self.all_faces_from_node(node=node_to_remove)
                if not faces_with_other_node:
                    """
                    Note to AJ: There could be faulty logic here because of the PyCharm warning but I am pretty sure it
                    is all good to go.
                    """
                    # This will remove the node AND any edges it belongs to.
                    self.design_graph.remove_node(node_to_remove)
                    # In here we also need to update the edge_divide_nodes list:
                    if node_to_remove in self.edge_divide_nodes.keys():
                        self.edge_divide_nodes.pop(node_to_remove, None)
                    other_keys_to_pop = []
                    for k, v in self.edge_divide_nodes.items():
                        if node_to_remove in v:
                            other_keys_to_pop.append(k)
                    for k in other_keys_to_pop:
                        self.edge_divide_nodes.pop(k, None)

        # Finally we reassign plotting faces and also remove this merge_list
        self.plotting_faces['other'] = all_plotting
        self.valid_to_merge.remove(faces)

    def find_triangles_containing_edge(self, edge: list, og_face: tuple) -> list:
        """
        This function will take in an edge and find and return a list of all triangles containing this edge.
        :param og_face: Original face we are dividing so we do not add this face to the list.
        :param edge: Edge we are dividing in a divide rule
        :return: List of triangles containing this edge that also need to be divided
        """
        all_triangles = self.plotting_faces['other']
        triangles_containing_edge = []
        for t in all_triangles:
            if edge[0] in t and edge[1] in t and t != og_face:
                triangles_containing_edge.append(t)

        return triangles_containing_edge

    def divide_similar_face(self, triangle_to_divide: tuple, new_node: int) -> list:
        """
        This function supports the "Divide_Face" Function and divides a similar face that is sharing an edge to generate
        a valid PLY file type.
        :param triangle_to_divide: Triangle face we are dividing
        :param new_node: New node that we need to connect up
        :return: List of new triangle vertices to add
        """
        # First we need to find which combination of nodes is NOT collinear. There should only be 1 pair:
        new_triangles = []
        for v1, v2 in itertools.combinations(list(triangle_to_divide), 2):
            if not self.colinearity_check(v1=self.xyz_from_node(v1), v2=self.xyz_from_node(v2),
                                          v3=self.xyz_from_node(new_node)):
                # Once this is found, we simply have the following first new triangle:
                new_triangles.append((v1, v2, new_node))
        # Also, to create a new edge, we need to find which vertex of the triangle_to_divide appears the most as that
        # will be the new edge with new_node:
        flat_list = [value for tpl in new_triangles for value in tpl]
        counter = Counter(flat_list)
        node_to_add_edge_to = counter.most_common(1)[0][0]
        # Add the edge:
        self.design_graph.add_edge(new_node, node_to_add_edge_to)
        return new_triangles

    def verify_coplanar_faces(self) -> None:
        """
        This function works alongside find_coplanar_faces. This function serves the purpose to ensure that we are NOT
        removing any binding region or NP faces from out design. It then verifies to ensure that the faces in these
        lists ARE in fact 100% in contact and coplanar.
        :return: Nothing, just updates a list.
        """
        new_merge_list = []
        invalid_faces_to_merge = self.plotting_faces['NP'] + self.plotting_faces['binding_region']
        for i in self.valid_to_merge:
            # If there is only 2 shared faces, we just need to make sure that both are not a binding_region AND that the
            # shared edge is NOT a terminal edge.
            new_merge = []
            for coplanar_face in i:
                if coplanar_face not in invalid_faces_to_merge:
                    new_merge.append(coplanar_face)
            # Now we go on to verify that if there is only 2 faces that 1) they still share an edge and 2) that edge is
            # not terminal:
            if len(new_merge) == 2:
                shared_edge = list(set(new_merge[0]) & set(new_merge[1]))
                if len(shared_edge) == 2:
                    if self.design_graph.nodes[shared_edge[0]]['terminal'] and \
                            self.design_graph.nodes[shared_edge[1]]['terminal']:
                        # If both nodes are terminal, then this is a terminal edge and we can NOT remove:
                        new_merge = []  # I just make this a blank value
                else:
                    # If the shared edge is NOT 2 values, then we can NOT merge the face as they are not coplanar faces
                    new_merge = []
            elif len(new_merge) > 2:
                check_multiple = []
                for a in new_merge:
                    touching_faces = self.find_coincident_faces(face=a)
                    touching_faces.append(a)
                    if not all(item in touching_faces for item in new_merge):
                        # If this is False that means all of the faces may not be touching so we don't try and merge:
                        new_merge = []

            # Now, after doing this, if the new_merge is at least 2 long, we append to new_merge_list as these can
            # be merge together:
            if len(new_merge) >= 2:
                new_merge_list.append(new_merge)

        # After this loop, we can then reset the value for valid_to_merge:
        self.valid_to_merge = new_merge_list

    def merge_two_faces_V2(self, face1: tuple, face2: tuple) -> bool:
        """
        This function merges to coplanar faces together. It is different than merge_n_faces because the logic is
        different due to PLY requirements for the filetype.
        :return:
        """
        # First we need to find the shared edge of the passed in faces
        shared_edge_1 = set(face1) & set(face2)
        shared_edge_2 = set()
        # Since in a divide_face rule we split two sets of faces, we actually need to find the "pair" to this shared
        # edge. To do this, we grab the coincident triangles to each face:
        coincident_face1 = self.find_coincident_faces(face=face1)
        coincident_face2 = self.find_coincident_faces(face=face2)

        # Of these coincident faces, we know that 2 of them will share an edge. We start by removing the face1 and face2
        # since we know that these are obviously already in contact:
        face3, face4 = None, None
        coincident_face1.remove(face2)
        coincident_face2.remove(face1)
        for f1 in coincident_face1:
            for f2 in coincident_face2:
                common_values = set(f1) & set(f2)
                if len(common_values) >= 2:
                    face3 = f1
                    face4 = f2
                    shared_edge_2 = common_values
                    break
        if face3 is None or face3 == face4:
            # In this case, we can NOT merge this face without undoing graphical work that may not be possible
            # geometrically through some combination of other moves, so we return False meaning no merge possible
            return False

        # If we get here, we know the 4 total faces that should be replaced. We know that a node must be removed which
        # will remove its edges. The node to remove is the highest mode between shared_edge_1 and shared_edge_2:
        counter = Counter(list(shared_edge_2) + list(shared_edge_1))
        most_common_value = counter.most_common(1)[0][0]
        new_T1 = list(set(list(face1) + list(face2)))
        new_T2 = list(set(list(face3) + list(face4)))

        # Remove the most_common_value node:
        self.design_graph.remove_node(most_common_value)
        new_T1.remove(most_common_value)
        new_T2.remove(most_common_value)
        # Also check to see if most_common_value is in the edge_divide_nodes dictionary and pop there if needed:
        if most_common_value in self.edge_divide_nodes.keys():
            self.edge_divide_nodes.pop(most_common_value, None)
        other_keys_to_pop = []
        for k, v in self.edge_divide_nodes.items():
            if most_common_value in v:
                other_keys_to_pop.append(k)
        for k in other_keys_to_pop:
            self.edge_divide_nodes.pop(k, None)

        all_plot = self.plotting_faces['other']

        for nT in [new_T1, new_T2]:
            for v1, v2 in itertools.combinations(nT, 2):
                if not self.design_graph.has_edge(v1, v2):
                    self.design_graph.add_edge(v1, v2)  # Add new triangle edges
            # We also add this new triangles:
            all_plot.append(tuple(nT))
            self.all_faces[tuple(nT)] = MeshFaces(v1=self.xyz_from_node(nT[0]), v2=self.xyz_from_node(nT[1]),
                                                  v3=self.xyz_from_node(nT[2]))


        # Removing old information:
        for i in [face1, face2, face3, face4]:
            all_plot.remove(i)
            self.all_faces.pop(i, None)

        # Re-assign the plotting_faces:
        self.plotting_faces['other'] = all_plot
        return True

    ####################################################################################################################
    ########################################### Production Rules #######################################################
    ####################################################################################################################
    def generate_start_shape(self) -> None:
        all_verts = np.array(self.get_all_verts())  # Grab all the vertices
        tri = Delaunay(all_verts)
        hull = tri.convex_hull

        # Extract edges from the convex hull
        faces = []
        unique_faces = []
        first_pass = True
        for simplex in hull:
            combs = list(itertools.combinations(simplex, 2))
            for c in combs:
                ## FINDING EDGES:
                edge = tuple(sorted(c))
                if edge not in self.design_graph.edges() and not self.nodes_on_same_face(n1=edge[0], n2=edge[1]):
                    self.design_graph.add_edge(edge[0], edge[1])

                ## FINDING FACES:
                opposite_vertex = [v for v in simplex if v not in c][0]
                faces.append(tuple(list(c) + [opposite_vertex]))
                if first_pass:
                    unique_faces.append(tuple(list(c) + [opposite_vertex]))
                    first_pass = False
            first_pass = True

        # Now we add the faces to our list, we do it this way BECAUSE we dont want to add the same face twice. First we
        # get a list of all the current faces we have:
        all_faces = []
        rect_faces = []
        exterior_predefined = []
        for i in self.all_faces.keys():
            if len(i) == 4:
                # If it is a rectangular face, we actually append every combination of 3 vertices so we do not re-add
                # an already existing face:
                for face_combo in itertools.combinations(list(i), 3):
                    rect_faces.append(sorted(list(face_combo)))
            else:
                all_faces.append(sorted(list(i)))

        for face in faces:
            if sorted(list(face)) not in all_faces and sorted(list(face)) not in rect_faces:
                self.plotting_faces['other'].append(face)
                self.all_faces[face] = MeshFaces(v1=self.xyz_from_node(node=face[0]),
                                                 v2=self.xyz_from_node(node=face[1]),
                                                 v3=self.xyz_from_node(node=face[2]), face_type='triangle')
                # We also append to make sure we don't read the same face:
                all_faces.append(sorted(list(face)))
            elif sorted(list(face)) in all_faces+rect_faces:
                # In this case, we want to actually store these nodal values as these are vertices used on the "outside"
                # of the Convex Hull:
                for i in sorted(list(face)):
                    if i not in exterior_predefined:
                        exterior_predefined.append(i)

        # Next we want to verify the unique faces and see if there are any edges (and therefore faces) to remove:
        og_faces = list(self.plotting_faces['NP'])
        keep_faces = []
        for uf in unique_faces:
            try:
                if all(item in self.og_NP_IDs for item in list(uf)):
                    # If this is true, we have a face from the NP that is on the Convex Hull, so we add this face. Note that
                    # due to triangulation and having square faces, we need to check to verify it is unique:
                    for np_face in list(self.plotting_faces['NP']):
                        if uf[0] in list(np_face) and uf[1] in list(np_face) and uf[2] in list(np_face) and np_face in \
                                og_faces:
                            og_faces.remove(np_face)
                            keep_faces.append(np_face)
            except AttributeError:
                pass

        # Now we need to remove the og_faces that are NOT a part of the Convex Hull:
        all_edges = []
        for f in all_faces:
            all_edges.append(sorted([f[0], f[1]]))
            all_edges.append(sorted([f[1], f[2]]))
            all_edges.append(sorted([f[0], f[2]]))

        removed_edges = []
        for i in og_faces:
            self.all_faces.pop(i, None)
            # After removing the face, we also need to verify the edges of the face whose pairings are below.
            NP_faces = [sorted([i[0], i[1]]), sorted([i[1], i[2]]), sorted([i[2], i[3]]), sorted([i[0], i[3]])]
            for f in NP_faces:
                if f not in all_edges and f not in removed_edges:
                    # If this edge from the removed face is NOT part of any edge of the new faces, we remove it from the
                    # design graph simply:
                    self.design_graph.remove_edge(f[0], f[1])
                    removed_edges.append([f[0], f[1]])

        # Finally, we re-assign the NP plotting faces:
        self.plotting_faces['NP'] = keep_faces

        # Next, we need to remove any un-needed nodal points which we can find via some lists:
        exterior_predefined = sorted(exterior_predefined)
        # Find the difference between the nodes used in the initial geoemtry from the convex hull generation
        dif = list(set(self.initial_nodes_pre_start) - set(exterior_predefined))
        # We simply remove the nodes found in the differences from the design graph:
        for i in dif:
            self.design_graph.remove_node(i)

        # Removed this, i do not think we need it but leaving here just in case i am stupid:
        """if self.routing_software == 'DAEDALUS':
            warnings.warn('NOTE: DAEDALUS and TALOS Edge Scaling is not yet implemented so generated structures may'
                          ' be much larger than expected.')"""

        # We also make a copy of the design_graph to use with the comparison:
        self.comparison_graph = deepcopy(self.design_graph)  # Make a copy

        # After assigning the comparison_graph, we also need to re-calculate the max number of basepairs:
        self.recalculate_max_number_bps()

        # Also, we need to calculate out the terminal edges, if any:
        terminal_faces = self.plotting_faces['NP'] + self.plotting_faces['binding_region']
        all_terminal_edges = []
        for f in terminal_faces:
            if len(f) == 4:
                terminal_edges = [(f[0], f[1]), (f[1], f[2]), (f[2], f[3]), (f[3], f[0])]
            else:
                terminal_edges = [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])]
            all_terminal_edges.extend(terminal_edges)

        for i in self.binding_regions:
            if i.binding_type == 'edge':
                p1, p2 = self.node_from_xyz(xyz=i.v1), self.node_from_xyz(xyz=i.v2)
                all_terminal_edges.append((p1, p2))
        self.terminal_edges = all_terminal_edges

    def recalculate_max_number_bps(self) -> None:
        """
        This function is called after the comparison graph is created. This will reassign the max number of basepairs
        so that we reject designs which exceed this number.
        :return: Nothing
        """
        # We go over all the edges of the comparison_graph and calculate these edge lengths depending on routing:
        total_length_in_nm = 0.
        for v1, v2 in self.comparison_graph.edges:
            length_in_nm = self.calculate_edge_length(n1=v1, n2=v2)
            total_length_in_nm += length_in_nm

        # We then divide dna rise to get an integer number of basepairs:
        numBasepairs = total_length_in_nm / self.dna_pitch_per_rise

        # Now if we are using vHelix, the "Ideal" software we are chasing is DAEDALUS (i.e. 2*nBaseapirs). If we are
        # routing with DAEDALUS, then the "Ideal" software is TALOS (i.e. 6*nbasepairs)
        if self.routing_software == 'vHelix':
            self.max_number_basepairs = 2 * int(numBasepairs)
        elif self.routing_software == 'DAEDALUS':
            self.max_number_basepairs = 6 * int(numBasepairs)

    def production_rules(self, rule_selection='') -> None:
        """
        This function is responsible for taking in the selected rule application and calling the proper function.
        :param rule_selection: String of function we are calling
        :return: Nothing
        """
        global design_constraint_failure  # Initialize this just in case a production rule does not apply
        if rule_selection == 'Divide Face':
            check = self.divide_face_rule()
        elif rule_selection == 'Merge Face':
            check = self.merge_face_rule()
        elif rule_selection == 'Extend Vertex':
            check = self.extend_vertex_rule()
        elif rule_selection == 'Retriangulate Face':
            check = self.retriangulate_face_rule()
        elif rule_selection == 'Normal Direction':
            check = self.move_vert_in_average_normal_dir_rule()
        else:
            raise Exception(f'Invalid rule passed in: {rule_selection}')

        # Now, if check is returned as False, we had an invalid rule applied so we need ot set design_constraint_failure
        if not check:
            design_constraint_failure = True

    def divide_face_rule(self) -> bool:
        """
        This function will take a face from the Convex Hull (as long as it is not a terminal face) and "divide it" into
        two smaller triangles. The "Check Edge Lengths" Design Constraint check will be very important regarding the
        validity / application of this rule.
        :return: If it returns True then all went well. If it returns False then this fails a design check!
        """
        all_verts = self.get_all_verts()
        # First we need to select a face that is not part of the binding regions or NP faces:
        face_verts_to_divide = choice(self.plotting_faces['other'])
        face_to_divide = self.all_faces[face_verts_to_divide]

        # Now we call the function from MeshFaces which will use an algorithm "divide" this face by using a bisection:
        newEdgePoint, newPoint, divided_edge, newT1, newT2 = face_to_divide.divide_face()

        # Now, if the newPoint we calculated IS ALREADY A POINT, then we have a non-planar, but touching (and symmetric)
        # face. This makes  graph modification easier as we just need to essentially add an edge (and not a new point)
        if np.any(np.all(newPoint == np.array(all_verts), axis=1)):
            # In this case, we simply first add a new edge connecting the newPoint to newEdgePoint
            n1, n2 = self.node_from_xyz(xyz=newEdgePoint), self.node_from_xyz(xyz=newPoint)
            self.design_graph.add_edge(n1, n2)

            # Next we remove the old faces:
            self.all_faces.pop(face_verts_to_divide, None)
            og_faces = self.plotting_faces['other']
            og_faces.remove(face_verts_to_divide)

            # And adding the new faces:
            for i in [newT1, newT2]:
                new_nodes = []
                for v in i:
                    # Loop over each vertex and find these nodes
                    new_nodes.append(self.node_from_xyz(xyz=v))
                # Adding in the new triangle to other and creating its mesh face:
                og_faces.append(tuple(new_nodes))
                self.all_faces[tuple(new_nodes)] = MeshFaces(v1=i[0], v2=i[1], v3=i[2], face_type='triangle')

            # Finally we reassign the value for plotting faces to our new og_faces list and return True meaning we
            # applied rule successfully so we don't do the set of code outside if.
            self.plotting_faces['other'] = og_faces
            return True

        # In the case the If is false, we have a unique, non-planar, non-symmetrical face we are dividing for the
        # first time so:

        # To begin modifying the graph, first we check the edge we are adding a point to to tell if it is terminal or
        # not:
        edge_verts = [self.node_from_xyz(divided_edge[0]), self.node_from_xyz(divided_edge[1])]
        terminal_node = False  # Start off by assuming the new node is NOT terminal
        t_count = 0
        for e in edge_verts:
            if self.design_graph.nodes[e]['terminal'] and 'binding_zone' not in self.design_graph.nodes[e]:
                t_count += 1
            if t_count == 2:
                # If both checked nodes are ONLY terminal points then the new node IS terminal
                terminal_node = True
        # Now we remove the edge in edge verts. We only need to remove this edge IF it is still an edge. Sometimes the
        # divide_face rule may get rid of an edge if faces share these.
        try:
            self.design_graph.remove_edge(edge_verts[0], edge_verts[1])
        except nx.NetworkXError:
            self.display_graph()
            raise Exception('There is an issue somewhere in the code logic. This should NOT normally raise...')
        # Add in the new vertex:
        self.add_new_node(x=newPoint[0], y=newPoint[1], z=newPoint[2], terminal=terminal_node)
        # We then add in the new edges. The add_new_node increments the node-count which is why we take -1 below:
        self.design_graph.add_edge(edge_verts[0], self.node_count-1)
        self.design_graph.add_edge(edge_verts[1], self.node_count-1)
        # Also add this to a master dictionary to use for comparisons later:
        self.edge_divide_nodes[self.node_count-1] = edge_verts
        # Adding in the new edge to create the two triangles:
        self.design_graph.add_edge(self.node_from_xyz(xyz=newEdgePoint), self.node_count-1)

        # Then, we must remove the face_to_divide and replace it with the two new triangles
        self.all_faces.pop(face_verts_to_divide, None)
        og_faces = self.plotting_faces['other']
        og_faces.remove(face_verts_to_divide)
        # Adding the new triangles:
        for i in [newT1, newT2]:
            new_nodes = []
            for v in i:
                # Loop over each vertex and find these nodes
                new_nodes.append(self.node_from_xyz(xyz=v))
            # Adding in the new triangle to other and creating its mesh face:
            og_faces.append(tuple(new_nodes))
            self.all_faces[tuple(new_nodes)] = MeshFaces(v1=i[0], v2=i[1], v3=i[2], face_type='triangle')

        # NEXT: Since the original edge_verts may be used in more triangles, we need to find any triangles containing
        # this edge:
        divide_these = self.find_triangles_containing_edge(edge=edge_verts, og_face=face_verts_to_divide)
        # For all of these triangles, we need to find which of the vertices in the triangle is NOT collinear with the
        # newest "divide point" so that we can create NEW triangles.
        for t in divide_these:
            # In this loop, we essentially call a function to divide this triangle into 2 using the new node
            new_triangles = self.divide_similar_face(triangle_to_divide=t, new_node=self.node_count-1)
            # Then we remove t from og_faces and add these new triangles to og_faces
            og_faces.remove(t)
            self.all_faces.pop(t, None)
            for new_t in new_triangles:
                og_faces.append(new_t)
                xyz = []
                for i in list(new_t):
                    xyz.append(self.xyz_from_node(i))
                self.all_faces[new_t] = MeshFaces(v1=xyz[0], v2=xyz[1], v3=xyz[2], face_type='triangle')

        # Lastly, reassign the plotting faces:
        self.plotting_faces['other'] = og_faces

        # If we get here, we return True meaning the rule applied sucessfully:
        return True

    def merge_face_rule(self) -> bool:
        """
        This function will combine coplanar faces to make one larger face. It will also be in charge of removing any
        "empty" nodes if they no longer belong.
        :return: If it returns True then all went well. If it returns False then this fails a design check!
        """
        # First, due to computational costs, we call the function to find the planar faces only here because we
        # do NOT want to constantly update this.
        self.find_coplanar_faces()
        self.verify_coplanar_faces()

        # Since we constantly update the coplanar faces, we look to see if there are any merge-able faces:
        if not self.valid_to_merge:
            # If there are no valid faces to merge, we just return False meaning INVALID rule choice!
            return False

        # Otherwise we randomly pick a pair of faces to merge. This will ALWAYS be a pair of two so:
        merge_choice = choice(self.valid_to_merge)
        if len(merge_choice) == 2:
            valid_merge = self.merge_two_faces_V2(face1=merge_choice[0], face2=merge_choice[1])
            if not valid_merge:
                return False
        else:
            self.merge_n_faces(faces=merge_choice)

        # If we get here without issue, we return True meaning it applied successfully
        return True


    def point_on_terminal_edge(self, selected_node: int) -> tuple[bool, tuple]:
        """
        This is used with the extend_vertex_rule and checks if a point is on a terminal edge. IF this is TRUE, then we
        can only move towards either of the terminal nodes
        """
        for e in self.terminal_edges:
            points_colinear = self.colinearity_check(v1=self.xyz_from_node(node=selected_node),
                                                     v2=self.xyz_from_node(node=e[0]),
                                                     v3=self.xyz_from_node(node=e[1]))
            if points_colinear:
                return True, e

        # If we loop through all and do not find a colinear point, we return False, (-1, -1) for a dummy value:
        return False, (-1, -1)


    def extend_vertex_rule(self) -> bool:
        """
        This rule will move a non-terminal node along the edge it is located.
        :return: False: Rule not applied due to no nonterminal node; True: Applied properly
        """
        # First we update the list of nonterminal nodes, and if there are not any we return False
        self.find_nonterminal_nodes()
        if not self.nonterminal_nodes:
            # If there are no nonterminal nodes, we return False meaning no valid rule application
            return False

        # Otherwise, we start by getting all faces referencing this node:
        selected_node = choice(self.nonterminal_nodes)
        all_faces = self.all_faces_from_node(node=selected_node)

        # Before we continue in this rule application, we need to see if the faces are coplanar or not as that will
        # change the direction we can extend to:
        ref_pt = self.get_3D_reference()
        #coplanar_faces = self.group_coplanar_faces(faces=all_faces, ref_pt=ref_pt)

        # If the selected node was in a previous edge division, we also want to verify other nodes that could be causing
        # issues:
        reassign = True
        vert_to_move_towards = None  # initializing to get rid of annoying warning from Pycharm.
        if selected_node in self.edge_divide_nodes.keys():
            # First we find what edges are there:
            nodes_on_line = self.find_nodes_on_line(node=selected_node)
            faces_on_edge = self.find_triangles_via_division(verts=nodes_on_line)
            cp = self.group_coplanar_faces(faces=(all_faces + faces_on_edge), ref_pt=ref_pt)
            # Now we just add this to the list of all_faces we are checking for.
            check_me = all_faces + faces_on_edge
            # Now if the lengths of coplanar faces and check_me are NOT equal, then we have a node that is on an edge that
            # we need to limit the direction it can move in:
            if len(cp) != len(check_me):
                # In this case, we must move in the direction of the nodes on the line, so we just simply choose one
                # randomly:
                vert_to_move_towards = choice(nodes_on_line)
                reassign = False

        if reassign:
            # Next we select a face at random (to establish a directionality) and find which vertex is the selected_node
            # and STORE THIS POINT (as we will need it later)
            change_face = list(choice(all_faces))
            change_face.remove(selected_node)  # Remove the selected node as a possibility
            # We also choose a random "node" to move towards on this face:
            vert_to_move_towards = choice(change_face)

        # Lastly, we actually re-assign the value of vert_to_move_towards if this nonterminal node is on a terminal edge
        on_term_edge, edge_of_point = self.point_on_terminal_edge(selected_node=selected_node)
        if on_term_edge:
            vert_to_move_towards = choice(edge_of_point)

        og_xyz_to_search_for = self.xyz_from_node(node=selected_node)
        xyz_moving_towards = self.xyz_from_node(node=vert_to_move_towards)

        # Find the encompassing vectors to determine which "direction" to move the vector along:
        u_v = (xyz_moving_towards - og_xyz_to_search_for) / np.linalg.norm((xyz_moving_towards - og_xyz_to_search_for))

        # Calculate the new nodal value for our selected_node and update the graph values and the face values:
        new_value = og_xyz_to_search_for + (self.extend_rule_distance * u_v)
        self.update_nodal_position(node=selected_node, new_xyz=new_value)
        self.update_face_values(node=selected_node, old_xyz=og_xyz_to_search_for, new_xyz=new_value)


        # If we get here without issue, we return True meaning it applied successfully
        return True

    def retriangulate_face_rule(self, override = None) -> bool:
        """
        This function choose a face in the mesh and refines it. This rule most likely has a small success rate as it
        will inherenetly create "small" faces that may not fit within design constraint but it is very useful for
        large faces / open designs
        :return: True: Rule applied successfully; False: Not applied
        """
        # First we need to select a face that is not part of the binding regions or NP faces:
        face_verts = choice(self.plotting_faces['other'])
        if override is not None:
            face_verts = override
        face_to_triangulate = self.all_faces[face_verts]
        # Finding the face centroid:
        centroid = face_to_triangulate.calculate_face_centroid()

        # Next we remove the original faces:
        self.remove_faces(face_verts=face_verts)

        # Then we simply add a new vertex here and create edges to the other face verts
        self.add_new_node(x=centroid[0], y=centroid[1], z=centroid[2], terminal=False)
        # NOTE: we add edges to node_count-1 because the add_new_node above will add 1 to the node_count.
        self.design_graph.add_edge(face_verts[0], self.node_count-1)
        self.design_graph.add_edge(face_verts[1], self.node_count-1)
        self.design_graph.add_edge(face_verts[2], self.node_count-1)

        # Finally we add the 3 new faces to all_faces and plotting_faces['other']:
        all_faces = self.plotting_faces['other']
        for v1, v2 in itertools.combinations(list(face_verts), 2):
            all_faces.append((self.node_count-1, v1, v2))
            v11 = self.xyz_from_node(node=self.node_count-1)
            v22 = self.xyz_from_node(node=v1)
            v33 = self.xyz_from_node(node=v2)
            self.all_faces[(self.node_count-1, v1, v2)] = MeshFaces(v1=v11, v2=v22, v3=v33)

        # If we get here without issue, we return True meaning it applied successfully
        return True

    def move_vert_in_average_normal_dir_rule(self) -> bool:
        """
        This rule will select a non-terminal node and move it in the average normal direction of the faces that the
        node belongs to. It will randomly choose which of the two average normal vectors to travel within since the
        cross products are all taken using a similar notation.
        :return:
        """
        # First we update the list of nonterminal nodes, and if there are not any we return False
        self.find_nonterminal_nodes()
        if not self.nonterminal_nodes:
            # If there are no nonterminal nodes, we return False meaning no valid rule application
            return False

        # Otherwise, we start by getting all faces referencing this node:
        while_count_check = 0
        selected_node = -1
        while self.nonterminal_nodes:
            selected_node = choice(self.nonterminal_nodes)
            point_on_term_edge, ed = self.point_on_terminal_edge(selected_node=selected_node)
            if not point_on_term_edge:
                break  # If the random choice is not on a terminal edge, we break and use this selected node.
            # We remove nodes from the selection criteria as we go since this list is constantly re-updated with each
            # rule selection
            self.nonterminal_nodes.remove(selected_node)
            while_count_check += 1
            if while_count_check > 100:
                warnings.warn('The move_vert_in_normal_direction rule is not able to find a nonterminal node that is '
                              'not on a terminal edge. There may be an issue in your final structure to check.')
                return False

        all_faces = self.all_faces_from_node(node=selected_node)
        # Before we continue in this rule application, we need to see if the faces are coplanar or not as that will
        # change the direction we can extend to:
        ref_pt = self.get_3D_reference()

        # To determine the direction we need to move in, we loop over all faces the node belong to and mark the
        # normal vector directions:
        sum_normal1, sum_normal2 = np.zeros(3), np.zeros(3)
        for face in all_faces:
            norm1, norm2 = self.all_faces[face].return_both_normal_vectors(reference_pt=ref_pt)
            sum_normal1 += norm1
            sum_normal2 += norm2

        # Now, we want the average normal direction to be a normal so:
        average_normal_1 = sum_normal1 / np.linalg.norm(sum_normal1)
        average_normal_2 = sum_normal2 / np.linalg.norm(sum_normal2)

        # We choose one of these two directions
        u_v = choice([average_normal_1, average_normal_2])

        # Now update the positions and things:
        og_xyz_to_search_for = self.xyz_from_node(node=selected_node)
        new_value = og_xyz_to_search_for + (self.extend_rule_distance * u_v)
        # Then we update the positions for this node to all lists:
        self.update_nodal_position(node=selected_node, new_xyz=new_value)
        self.update_face_values(node=selected_node, old_xyz=og_xyz_to_search_for, new_xyz=new_value)
        # And return True meaning the rule applied properly:
        return True


    ####################################################################################################################
    ######################################## Design Constraint Checking ################################################
    ####################################################################################################################
    def master_design_constraint_checker(self) -> bool:
        """
        This function is called every internal loop of Shape Annealing to verify that the moves we are making are
        leading to a valid design based on my definitions OR software constraints.

        :return: True: VALID design, False: INVALID Design
        """
        c1 = self.outside_design_space()
        c2 = self.vert_in_NP_zone()
        c3 = self.edge_in_NP_zone()
        c4 = self.invalid_edge_lengths()
        c5 = self.invalid_production_rule()
        c6 = self.used_too_much_scaffold()
        c7 = self.invalid_face_angles()
        c8 = self.too_many_edges_at_a_node()
        c9 = self.invalid_mesh()
        if self.allow_self_intersecting_geometry:
            c10 = self.intersecting_faces()  # We always assume the faces are not intersecting here!
        else:
            #c10 = self.intersecting_faces()
            c10 = self.find_intersecting_edges()

        # Now we just see if any of these checks are True meaning there was an invalid move taken leading to an
        # invalid design:
        if self.check_scaffold_length_constraint:
            all_checks = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10]
        else:
            # If the above is turned off, we do not check the used_too_much_scaffold check.
            all_checks = [c1, c2, c3, c4, c5, c7, c8, c9, c10]
        if any(all_checks):
            # If any of these values are True we have an INVALID design and return False meaning the design is invalid:
            return False
        else:
            # Otherwise we return True meaning GOOD (or non-design constraint breaking) design.
            return True

    def outside_design_space(self) -> bool:
        """
        This function checks to see if anything in the generative design space is "outside" the limits of the input
        unit cell parameters from the user.
        :return: True: Design is INVALID!!!!, False: Design Constraint Check passed
        """
        # First we get all vertices:
        all_verts = np.array(self.get_all_verts())
        all_X = all_verts[:, 0]
        all_Y = all_verts[:, 1]
        all_Z = all_verts[:, 2]

        # Now we just run a few if statement checks:
        result_X = np.any((all_X < 0.) | (all_X > self.maxX))
        result_Y = np.any((all_Y < 0.) | (all_Y > self.maxY))
        result_Z = np.any((all_Z < 0.) | (all_Z > self.maxZ))

        # If any of the above are True, then we return True:
        if any([result_X, result_Y, result_Z]):
            return True

        # Otherwise we return False meaning the rule applied properly:
        return False

    def vert_in_NP_zone(self) -> bool:
        """
        This function checks to verify that we are not entering areas that are "Do Not Enter" regions for a user
        defined porosity reason, for example.
        :return: True: Design is INVALID, False: Design is VALID
        """
        all_verts = self.get_all_verts()
        # We must loop over all NP regions (and eventually the NP Regions)
        for r in self.NP_Regions:
            # Now we loop over all the vertices to see if they are inside this region.
            for v in all_verts:
                check = r.point_inside_space(v)
                if check:
                    # If we are ever inside the NP Zone space, we just return True meaning INVALID DESIGN
                    return True

        # If we get here, we return False meaning the rule applied properly:
        return False

    def edge_in_NP_zone(self) -> bool:
        """
        This is similar to above, but we are now just checking if any edge intersects one of the faces of an NP Zone
        :return: True: Design is INVALID, False: Design is VALID
        """
        # We must loop over all NP regions (and eventually the NP Regions)
        checked_edges = []
        for r in self.NP_Regions:
            # Now we loop over all the vertices to see if they are inside this region.
            for v1, v2 in self.design_graph.edges():
                if (v1, v2) not in checked_edges:
                    # First we get the XYZ coordinates of each point:
                    p1, p2 = self.xyz_from_node(node=v1), self.xyz_from_node(node=v2)
                    checked_edges.extend([(v1, v2), (v2, v1)])  # Append both edge combs to ensure we don't double check
                    # Now we just pass to see:
                    check = r.edge_intersecting_face(p1=p1, p2=p2)
                    if check:
                        # If we are ever inside the NP Zone space, we just return True meaning INVALID DESIGN
                        return True

        # If we get here, we return False meaning the rule applied properly:
        return False

    def invalid_edge_lengths(self) -> bool:
        """
        This function is there to ensure none of the edges we are generating are below (or above) a certain length.
        :return: True: Design is INVALID, False: Design is VALID
        """
        all_lengths = []
        checked_edges = []
        for v1, v2 in self.design_graph.edges():
            if (v1, v2) not in checked_edges:
                # First we get the XYZ coordinates of each point:
                p1, p2 = self.xyz_from_node(node=v1), self.xyz_from_node(node=v2)
                dist = self.calculate_3D_distance(p1=p1, p2=p2)
                all_lengths.append(dist)
                checked_edges.extend([(v1, v2), (v2, v1)])
        # Reassign to a NumPy array for simple checking:
        all_lengths = np.array(all_lengths)
        # If any edge is below the min edge length threshold or above the max edge length threshold we return True
        ## First we need to convert from basepairs to nanometers for checking:
        if self.routing_software == 'DAEDALUS':
            minEdge = round(self.dna_pitch_per_rise * self.min_edge_length)
            maxEdge = round(self.dna_pitch_per_rise * self.max_edge_length)
        else:
            minEdge = self.min_edge_length
            maxEdge = self.max_edge_length
        if np.any((all_lengths < minEdge) | (all_lengths > maxEdge)):
            return True

        # If we get here, we return False meaning the rule applied properly:
        return False

    def used_too_much_scaffold(self) -> bool:
        """
        A simple design constraint check that rejects designs that are estimated to use more than the user defined input
        number of basepairs
        :return: True: Design is INVALID, False: Design is VALID
        """
        numBPs = self.calculate_total_edge_lengths()
        if numBPs > self.max_number_basepairs:
            return True
        else:
            return False

    def invalid_production_rule(self) -> bool:
        """
        This function checks to see if any of the production rules were invalidly applied. For example, if we do not
        have any faces to merge but Shape Annealing selected a "Merge Face" rule, then we need to catch this.
        :return: True: Design is INVALID, False: Design is VALID
        """
        global design_constraint_failure
        if design_constraint_failure:
            # If this is true we failed to apply a rule and therefore reject this design iteration (and reset the flag)
            design_constraint_failure = False
            return True
        else:
            # Otherwise we just return False meaning this design constraint check passed and there is no issues.
            return False

    def invalid_face_angles(self) -> bool:
        """
        This function looks at all the mesh faces (that are triangular) and looks for the ones < the threshold
        :return: True: Design is INVALID, False: Design is VALID
        """
        mesh_file = self.invalid_mesh(return_mesh=True)

        for angles in mesh_file.face_angles:
            face_angles = []
            for i in angles:
                face_angles.append(np.rad2deg(i))
            if any(value < self.min_face_angle for value in face_angles):
                return True

        """for k, v in self.all_faces.items():
            # We only check the triangular faces since it's all we care about:
            if len(k) == 3:
                # First we recalculate the face angles:
                v.calculate_face_angles()
                all_angles = v.face_angles
                if any(value < self.min_face_angle for value in all_angles):
                    return True"""

        # If we get to the end, that means all face angles ARE valid so we return False:
        return False

    def too_many_edges_at_a_node(self) -> bool:
        """
        There is a limit to how many edges can be at any given node. Having too many edges at one node can lead to an
        increase in the average RMSF of a given structure. This uses a user defined input for "max_number_edges_at_node"
        :return: True: Design is INVALID, False: Design is VALID
        """
        for n in self.design_graph.nodes:
            temp = self.design_graph.degree(n)
            if self.design_graph.degree(n) > self.max_number_edges_at_node:
                # If any node has a degree (i.e. # of edges) >= the set value, we return False for invalid design
                return True
        # If we check all nodes and never return True, then we return False meaning our design is valid.
        return False

    def invalid_mesh(self, return_mesh: bool = False) -> Union[bool, trimesh.Trimesh]:
        """
        This function gurantees that a mesh is "valid" and can be validates
        :return: True: Design is INVALID, False: Design is VALID
        """
        # Gather and Reformat Vertices and Face information needed:
        vertices = np.array(self.get_all_verts())
        faces = np.array(list(self.all_faces.keys()), dtype=object)

        # Renumber face values:
        vertex_mapping = self.renumber_vertices(faces=faces)
        reference_point = np.mean(vertices, axis=0)

        all_faces = []
        for row in faces:
            # We first use the mapping dictionary to ensure we are accessing the correct index positions:
            new_row = []
            for ind in list(row):
                new_row.append(vertex_mapping[ind])
            # Before writing the face, we need to process it to verify that it is "correct":
            face_vertices = vertices[np.array(new_row)]
            face_normal = self.calculate_face_normal(face_vertices)

            if not self.is_face_clockwise(face_vertices, face_normal, reference_point):
                # If the face is clockwise, we re-orient properly depending on if rectangular / triangular face:
                if len(new_row) == 4:
                    new_row1 = [new_row[0], new_row[3], new_row[2]]
                    new_row2 = [new_row[0], new_row[2], new_row[1]]
                    all_faces.append(new_row1)
                    all_faces.append(new_row2)
                else:
                    new_row_1 = [new_row[0], new_row[2], new_row[1]]
                    all_faces.append(new_row_1)
            # If we didn't change anything:
            else:
                if len(new_row) == 4:
                    new_row1 = [new_row[0], new_row[1], new_row[2]]
                    new_row2 = [new_row[0], new_row[2], new_row[3]]
                    all_faces.append(new_row1)
                    all_faces.append(new_row2)
                else:
                    new_row_1 = [new_row[0], new_row[1], new_row[2]]
                    all_faces.append(new_row_1)

        all_verts = []
        for v in vertices:
            all_verts.append(list(v))
        my_mesh = trimesh.Trimesh(vertices=all_verts, faces=all_faces)
        # Now fix the normals:
        my_mesh.fix_normals()
        if return_mesh:
            return my_mesh
        else:
            if my_mesh.is_volume:
                return False  # If the mesh is a volume, then the design is VALID and we return False per our notation
            else:
                return True  # If the mesh is not a volume, we return True signifying "INVALID DESIGN!"

    def intersecting_faces(self) -> bool:
        """
        We should not have overlapping faces in a mesh file... maybe
        :return: True: Design is INVALID, False: Design is VALID
        """
        mesh = self.invalid_mesh(return_mesh=True)
        if mesh.is_watertight: # If the mesh is watertight we assume it is good to go and return False
            return False
        return True  # If the mesh is not watertight, we return True meaning invalid design

    def find_intersecting_edges(self) -> bool:
        """
        This function instead considers the design as an invalid polyhedra where faces may intersect one another. This
        specific function will consider edges has cylindrical bounding boxes that may not collide with one another
        through a collision detection algorithm.

        :return: True: Design is INVALID, False: Design is VALID
        """
        # To start, we find all pairs of 2 edges that are not connected in the graph.
        all_edges = list(self.design_graph.edges())
        for edge1 in self.design_graph.edges():
            edge_list_copy = deepcopy(all_edges)
            filtered_list = [tup for tup in edge_list_copy if edge1[0] not in tup]
            filtered_list = [tup for tup in filtered_list if edge1[1] not in tup]

            # Next, we check through all pairs of two un-connected edges with cylindrical box approximations to check if
            # there is any collisions in the converted DNA model
            num_points = 10  # This is for refine-ness (too many points may be overly restrictive)
            for edge2 in filtered_list:
                P1, P2 = self.xyz_from_node(node=edge1[0]), self.xyz_from_node(node=edge1[1])
                P3, P4 = self.xyz_from_node(node=edge2[0]), self.xyz_from_node(node=edge2[1])
                axis1 = P2 - P1
                axis2 = P4 - P3
                # Calculate the step size for distributing points along the axis
                t_values = np.linspace(0, 1, num_points)
                # Calculate the points along axis1 and axis2:
                points_on_axis1 = P1 + t_values[:, np.newaxis] * axis1
                points_on_axis2 = P3 + t_values[:, np.newaxis] * axis2

                # Calculate distances from all points on axis1 to axis2
                distances = np.linalg.norm(points_on_axis1 - points_on_axis2, axis=1)

                # Check if any distance is less than 2 * R
                if np.any(distances < self.dna_diameter):
                    # If at any point we find a collision, we return True meaning the design is invalid and stop
                    # the search
                    return True  # Colliding
        # If all pairs are checked and not colliding, we return False meaning the design is good to go.
        return False


    ####################################################################################################################
    ########################################## Cost Function Functions #################################################
    ####################################################################################################################
    def calculate_total_edge_lengths(self) -> int:
        """
        This function returns the (estimated) number of nucleotides that DAEDALUS would use to route this design.
        :return: Int value in # of BPs used in our analysis
        """
        def calculate_graph_edge_lengths(design_graph, atrail=False):
            checked_edges = []
            total_length_in_nm = 0.
            for v1, v2 in design_graph.edges():
                if not atrail:
                    if (v1, v2) not in checked_edges and (v2, v1) not in checked_edges:
                        length_in_nm = self.calculate_edge_length(n1=v1, n2=v2)
                        total_length_in_nm += length_in_nm
                        checked_edges.append((v1, v2))
                elif atrail:
                    length_in_nm = self.calculate_edge_length(n1=v1, n2=v2)
                    total_length_in_nm += length_in_nm
                    checked_edges.append((v1, v2))
            # We then divide dna rise to get an integer number of basepairs:
            numBasepairs = total_length_in_nm / self.dna_pitch_per_rise
            return numBasepairs

        if self.routing_software == 'vHelix':
            # When using vHelix, there is a relatively computationally expensive A-Trail calculation. To get around this
            # I will use NetworkX to calculate a Eulerian circuit:
            # First check if Eulerian, and if so, we just calculate the nBasepairs similar to DAEDALUS since every edge
            # is traversed once in this case:
            graph_copy = deepcopy(self.design_graph)
            if nx.is_eulerian(graph_copy):
                nBasepairs = calculate_graph_edge_lengths(graph_copy)
                return int(nBasepairs)
            elif not nx.is_eulerian(graph_copy) and nx.is_connected(graph_copy):
                    # We also only check when the graph is connected. This is like 99% of the time except...
                    new_graph = nx.eulerize(graph_copy)
                    nBasepairs = calculate_graph_edge_lengths(new_graph, atrail=True)
                    return int(nBasepairs)
            elif not nx.is_connected(graph_copy):
                # At the very first pass and check of graph creation, there will be 0 edges. We just essentially
                # return a value of 0 and ignore this because no data in the first pass is stored (which is the only
                # time this else is activated).
                return 0


        elif self.routing_software == 'DAEDALUS':
            # Although DAEDALUS will make some "changes" to the design, this is the best way to approximate how much
            # scaffold is used without running DAEDALUS at every design iteration:
            graph_copy = deepcopy(self.design_graph)
            nBasepairs = calculate_graph_edge_lengths(graph_copy)
            return 2 * int(nBasepairs)

        elif self.routing_software == 'TALOS':
            graph_copy = deepcopy(self.design_graph)
            nBasepairs = calculate_graph_edge_lengths(graph_copy)
            return 6 * int(nBasepairs)
        else:
            raise Exception('Invalid routing software.')

    def cross_sectional_area_approximation(self) -> float:
        """
        Depending on the routing software being used there will be a different cross-sectional area!
        :return:
        """
        if self.routing_software == 'DAEDALUS':
            A = (2 * self.dna_diameter) * self.dna_diameter  # Approximate DX bundle as rectangular area
        elif self.routing_software == 'vHelix':
            A = (np.pi / 4) * self.dna_diameter**2
        else:
            raise Exception('Invalid routing software passed in.')
        return A

    def calculate_cost_function(self) -> float:
        """
        This function will be used to calculate the cost function values for an iteration of the design generation
        :return: A floating point value of the objective function valuation with no reward or penalty
        """

        # For the CURRENT DESIGN GRAPH:
        curNum = 0.  # Current numerator value
        graph_copy = deepcopy(self.design_graph)
        ## If our routing_software is DAEDALUS, this is an easier calculation to estimate, simply as:
        if self.routing_software == 'DAEDALUS':
            A = self.cross_sectional_area_approximation()
            for u, v in self.design_graph.edges:
                p1, p2 = self.xyz_from_node(node=u), self.xyz_from_node(node=v)
                L = self.calculate_3D_distance(p1=p1, p2=p2)
                dL = L - self.persistance_length
                curNum += (A * (-1 * dL**3) / self.persistance_length)
        ## If our routing software is vHelix, this is a harder calculation to estimate since some edges may have 2
        ## helices per edges and others may have one. We use a NetworkX eulerize function to quickly estimate these
        ## values.
        elif self.routing_software == 'vHelix':
            A1 = self.cross_sectional_area_approximation()  # Area used if there is a 1 helix edge
            A2 = (2 * self.dna_diameter) * self.dna_diameter  # Area used if there is a 2 helix edge:
            ## If the graph is simply eulerian, we can do a similar calculation to the DAEDALUS calc above:
            if nx.is_eulerian(graph_copy):
                for u, v in self.design_graph.edges:
                    p1, p2 = self.xyz_from_node(node=u), self.xyz_from_node(node=v)
                    L = self.calculate_3D_distance(p1=p1, p2=p2)
                    dL = L - self.persistance_length
                    # In the case it is Eulerian then each edge is a 1-helix bundle so we use A1:
                    curNum += (A1 * (-1 * dL ** 3) / self.persistance_length)

            ## Otherwise, we need to find the edges that occur twice (for the A-trail routing approximation) to find
            ## which Area values to use in the objective function calculation:
            elif not nx.is_eulerian(graph_copy) and nx.is_connected(graph_copy):
                # We also only check when the graph is connected. This is like 99% of the time except...
                new_graph = nx.eulerize(graph_copy)
                # Now we need to calculate which edges in the multigraph are used twice:
                checked_edges = []
                for v1, v2 in new_graph.edges():
                    p1, p2 = self.xyz_from_node(node=v1), self.xyz_from_node(node=v2)
                    L = self.calculate_3D_distance(p1=p1, p2=p2)
                    dL = L - self.persistance_length
                    if (v1, v2) not in checked_edges:
                        curNum += (A1 * (-1 * dL ** 3) / self.persistance_length)
                        checked_edges.append((v1, v2))
                    else:
                        # If we are looping over the same edge twice, we first remove the initial contribution from the
                        # 1HB from curNum, then add on the curNum for a 2HB on this edge:
                        curNum -= (A1 * (-1 * dL ** 3) / self.persistance_length)
                        curNum += (A2 * (-1 * dL ** 3) / self.persistance_length)

            else:
                raise Exception('This should NOT ever happen, if it does need to debug AJ!')

        else:
            raise Exception('Invalid routing software, no idea how that happens here!')

        # Now, depending on if we are using vHelix or DAEDALUS, the denominator term changes a bit. This is because the
        # goal of vHelix is to mimic a design that is stable as a DAEDALUS design of the start shape. However, if we are
        # using DAEDALUS, then we want to maek a design as stable as the TALOS design:
        if self.routing_software == 'vHelix':
            A_star = (2 * self.dna_diameter) * self.dna_diameter  # If using vHelix then A* is DAEDALUS calculation
        elif self.routing_software == 'DAEDALUS':
            A_star = (np.pi / 4) * (6**2)  # For a 6HB the approximate circular diameter is 6nm
        else:
            raise Exception('Invalid routing_software choice!')
        # For the COMPARISON DESIGN GRAPH:
        curDen = 0.
        for u, v in self.comparison_graph.edges:
            p1, p2 = self.xyz_from_node(node=u), self.xyz_from_node(node=v)
            L = self.calculate_3D_distance(p1=p1, p2=p2)
            dL = L - self.persistance_length
            curDen += (A_star * (-1 * dL**3) / self.persistance_length)

        return curNum / curDen

    def get_all_edge_lengths(self) -> np.array:
        """
        This function will return a numpy array of all the current edge lengths in nanometers
        :return:
        """
        all_lengths = []
        checked_edges = []
        for v1, v2 in self.design_graph.edges():
            if (v1, v2) not in checked_edges:
                # First we get the XYZ coordinates of each point:
                p1, p2 = self.xyz_from_node(node=v1), self.xyz_from_node(node=v2)
                dist = self.calculate_3D_distance(p1=p1, p2=p2)
                all_lengths.append(dist)
                checked_edges.extend([(v1, v2), (v2, v1)])
        # Reassign to a NumPy array for simple checking:
        all_lengths = np.array(all_lengths)
        return all_lengths

    def print_shortest_edge_length(self) -> int:
        """
        Simply prints to the termimnal what the shortest edge length is
        :return:
        """
        global numDecimals
        shortest_edge = 1e6  # Initialize
        n1, n2 = -1, -1
        for u, v in self.design_graph.edges:
            dist = self.calculate_3D_distance(p1=self.xyz_from_node(node=u), p2=self.xyz_from_node(node=v))
            dist = round(dist, numDecimals)
            if dist < shortest_edge:
                n1 = u
                n2 = v
                shortest_edge = dist
        """print(f'The shortest edge is between node {n1} and {n2} with a distance of {shortest_edge} nm'
              f' or {int(shortest_edge / 0.34)} basepairs.')"""
        return int(shortest_edge / 0.34)

    def calculate_all_surface_areas(self) -> tuple[int, list]:
        """
        This function will calculate all of the surface areas of the current design to use in an objective function.
        It will also return the number of faces.
        :return:
        """
        faces_to_check = self.plotting_faces['other']  # We don't calculate the terminal geometry faces
        total_face_areas = []
        #face_data = []  # Will be used to track the areas of all faces.
        for f in faces_to_check:
            # Extract the area and add to the total:
            face = self.all_faces[f]
            face_area = face.recalculate_face_area()
            total_face_areas.append(face_area)
            #total_face_area.append(face.area)
            #face_data.append(face.area)
        #std_dev = float(np.std(np.array(face_data)))
        #return len(faces_to_check), total_face_area, std_dev
        return len(faces_to_check), total_face_areas

    def number_of_edges_per_node(self) -> float:
        totalN = 0
        numConnnections = 0
        for n in self.design_graph.nodes:
            numConnnections += self.design_graph.degree(n)
            totalN += 1
        return numConnnections/totalN

    def minimize_RBD_forces(self) -> tuple[float, float, int, int]:
        """
        This function is used in the cost function to minimize the total RBD forces in the overall system.
        Note to AJ: The weight factors / spring constants and things will be completed in the Shape_Anealing script
                    since these constants are largely used to "scale" our objective function value.
        :return:
        """
        def inertia(n):
            # The 2nd moment of inertia for an N-helix bundle is calculated as:
            temp = ((2 / np.pi**2) * n**3 + n) * (np.pi / 4) * (self.dna_diameter/2)**4
            return temp
        ### FIRST: Contributions from the overall edge lengths of the faces in the design as a buckling beam:
        F_buck = 0.
        N_e = len(list(self.comparison_graph.edges()))  # Initialize
        track_repulsion_values = []  # Used in tracking later
        if self.routing_software == 'DAEDALUS':  # If using daedalus, calculating the total buckling force is easier:
            N_e = len(list(self.design_graph.edges()))  # Number of edges
            for u, v in self.design_graph.edges:
                p1, p2 = self.xyz_from_node(node=u), self.xyz_from_node(node=v)
                L = self.calculate_3D_distance(p1=p1, p2=p2)
                # Assumes both ends of DNA are pinned, n=2 for 2HB in DAEDALUS:
                F_buck += (np.pi**2 * self.E_dsdna * inertia(2) / (L**2))

        ### SECOND: Contributions from the repulsion forces between the edges (which are "groups" in essence.
        # We want to loop over all edges exactly once so we don't double count any contributions:
        graph_copy = deepcopy(self.design_graph)
        F_rep = 0.  # Repulsion force
        for e1, e2 in itertools.combinations(list(self.design_graph.edges()), 2):
            # We first need the centroid of each edge pair, so starting with e1:
            center_e1 = np.mean([self.xyz_from_node(node=e1[0]), self.xyz_from_node(node=e1[1])], axis=0)
            center_e2 = np.mean([self.xyz_from_node(node=e2[0]), self.xyz_from_node(node=e2[1])], axis=0)
            # Now, we want to find the distance between these centroids for d:
            d = self.calculate_3D_distance(p1=center_e1, p2=center_e2)

            # Now, depending on the edge and routing software, the "radius" values for the groupings depends.
            if self.routing_software == 'DAEDALUS':
                # In this case, the radius for each bundle is CONSTANT, and thus is 2 helices * D_dna / 2 (for radius):
                r_a, r_b = self.dna_diameter, self.dna_diameter
            elif self.routing_software == 'vHelix':
                # This case is much trickier and non-trivial due to an A-trail. If the graph is Eulerian, it is easy:
                if nx.is_eulerian(graph_copy):
                    N_e = len(list(self.design_graph.edges()))
                    # If Eulerian, then all the edges are 1HB and therefore the radius of the group is 1/2 the diameter
                    # of a single helix of DNA
                    r_a, r_b = (self.dna_diameter / 2), (self.dna_diameter / 2)
                    # In this case, the buckling force is also easily calculated:
                    p1, p2 = self.xyz_from_node(node=e1[0]), self.xyz_from_node(node=e1[1])
                    p3, p4 = self.xyz_from_node(node=e2[0]), self.xyz_from_node(node=e2[1])
                    L1 = self.calculate_3D_distance(p1=p1, p2=p2)
                    L2 = self.calculate_3D_distance(p1=p3, p2=p4)
                    # Assumes E=1 (constant anyways), both ends of DNA are pinned, n=2 for 2HB in DAEDALUS:
                    F_buck += (np.pi ** 2 * inertia(1) * self.E_dsdna / (L1 ** 2))
                    F_buck += (np.pi ** 2 * inertia(1) * self.E_dsdna / (L2 ** 2))
                else:
                        # Otherwise, we need to create an A-trail graph and find which edges appear twice, and if that is
                        # the edge we are searching, then one (or both) of our radius values changes:
                        new_graph = nx.eulerize(graph_copy)
                        N_e = len(list(new_graph.edges()))
                        new_graph_edges = list(new_graph.edges())
                        e1_count = new_graph_edges.count(e1)
                        e2_count = new_graph_edges.count(e2)
                        # Here, we also need to calculate buckling force which is a bit more wonky:
                        if e1_count == 1:  # If an edge appears once, then it is a single HB:
                            r_a = self.dna_diameter / 2
                            # If e1 is a single helix bundle we add these contributions to buckling force:
                            p1, p2 = self.xyz_from_node(node=e1[0]), self.xyz_from_node(node=e1[1])
                            L1 = self.calculate_3D_distance(p1=p1, p2=p2)
                            F_buck += (np.pi ** 2 * inertia(1) * self.E_dsdna / (L1 ** 2))
                        else:  # Otherwise it is a 2HB edge:
                            r_a = self.dna_diameter
                            # If e1 is a 2HB, then we change the inertia term but rest stays the same:
                            p1, p2 = self.xyz_from_node(node=e1[0]), self.xyz_from_node(node=e1[1])
                            L1 = self.calculate_3D_distance(p1=p1, p2=p2)
                            F_buck += (np.pi ** 2 * inertia(2) * self.E_dsdna / (L1 ** 2))
                        if e2_count == 1:
                            r_b = self.dna_diameter / 2
                            # Similar to above just using e2 instead:
                            p1, p2 = self.xyz_from_node(node=e2[0]), self.xyz_from_node(node=e2[1])
                            L1 = self.calculate_3D_distance(p1=p1, p2=p2)
                            F_buck += (np.pi ** 2 * inertia(1) * self.E_dsdna / (L1 ** 2))
                        else:
                            r_b = self.dna_diameter
                            p1, p2 = self.xyz_from_node(node=e2[0]), self.xyz_from_node(node=e2[1])
                            L1 = self.calculate_3D_distance(p1=p1, p2=p2)
                            F_buck += (np.pi ** 2 * inertia(2) * self.E_dsdna / (L1 ** 2))
            else:
                raise Exception('Invalid routing software.')
            # Now, regardless of routing software, the repulsion is simply:
            repulsion_value = (self.repulsion_distance - (d / (r_a + r_b)))
            track_repulsion_values.append(repulsion_value)

            # If the repulsion value is negative or zero, we use a small exponential term to show there is a minimal
            # amount of repulsion between these two edges, Otherwise, we append the repulsion value:
            if repulsion_value <= 0:
                new_r = (d / (r_a + r_b)) - self.repulsion_distance # Get the inverse value of repulsion_value
                scale_exponentially = np.exp(-1 * new_r)  # Scaled to be between 0 and 1 since the min value of the
                # new_r will be 0
                # Now using the scaled term, we add a VERY SMALL repulsion (0.01 * scale) to increment repulsion:
                F_rep += 0.01 * scale_exponentially  # We add a very MINIMAL (but non-0) amount of repulsion to help
                # guide the space search.
            else:
                F_rep += (self.C_repulsion * repulsion_value)  # C_repulsion makes sure it scales!

        ### THIRD AND FINALLY: The last piece we need is the term based on number of (etimated) nucleotides in design i
        #all_faces = self.plotting_faces['other']
        numBasepairs = self.calculate_total_edge_lengths()
        return F_buck, F_rep, numBasepairs, N_e


    ####################################################################################################################
    ########################################## Final Graph Processing ##################################################
    ####################################################################################################################
    def renumber_vertices(self, faces: np.array) -> dict:
        """
        This function will take in the list of vertices and faces and re-number the faces to verify the index positions
        :return:
        """
        # First we need to get all the unique face indices:
        face_vertex_mapping = {}
        unique_faces = []
        for i in faces:
            for v in list(i):
                if v not in unique_faces:
                    unique_faces.append(v)
        # Next we sort this list and start "counting" through the list to create a mapping via a dictionary:
        unique_faces = sorted(unique_faces)
        vertex_counter = 0
        for i in unique_faces:
            face_vertex_mapping[i] = vertex_counter
            vertex_counter += 1
        return face_vertex_mapping

    def calculate_face_normal(self, vertices: np.array) -> np.array:
        """
        This function is used to calculate the cross-product of a face so we can ensure we are writing the normal dir
        :param vertices: Vertex Points used
        :return:
        """
        global numDecimals
        edge1 = vertices[1] - vertices[0]
        edge2 = vertices[-1] - vertices[0]
        return np.cross(edge1, edge2)

    def is_face_clockwise(self, vertices: np.array, face_normal: np.array, reference_point: np.array) -> bool:
        """
        This function will take in a faces vertices and the calculated normal direction of the face and determine if
        the vertices are in a clockwise or counter clockwise rotation
        :param vertices: List of vertices xyz coordinates for a face
        :param face_normal: Cross product of two connected edges of the face
        :param reference_point: A standardized location to take dot products from
        :return: True: CCW vertices listed; False: Clockwise and need to re-orient the faces
        """
        centroid_to_reference = reference_point - np.mean(vertices, axis=0)
        dot_product = np.dot(face_normal, centroid_to_reference)
        return dot_product < 0

    def write_ply_file(self, filepath: str = None, filepath2: str = None, filename: str = 'Generated_Design'):
        """
        This function will write the ply-compliant file of the various faces in our 3D file.

        NOTE: If you are reading this, the write_p`ly_file is an insane amount of spaghetti code that I just didn't
              bother fixing up. Ignore this all, it should work!!!

        :param filepath: Where we are saving the PLY file to
        :param filepath2: Where we are saving the .aj file for graph comparison metrics
        :param filename: Name of PLY file WITHOUT the file extension
        :return: Nothing, just writes a plyfile
        """
        if filepath is None:
            filename2 = filename
            # If no filepath specified, we just write to the local SavedPlyFiles default directory:
            saveDir = os.getcwd() + '/SavedPlyFiles/'
            if not os.path.exists(saveDir):
                # Check to create a filepath if needed.
                os.makedirs(saveDir)
            if '.ply' not in filename:
                filename += '.ply'
                filename2 += '.aj'

            writeDir = saveDir + filename
            writeDir2 = saveDir + filename2
        else:
            writeDir = filepath
            writeDir2 = filepath2

        # Gather and Reformat Vertices and Face information needed:
        vertices = np.array(self.get_all_verts())
        faces = np.array(list(self.all_faces.keys()), dtype=object)

        # Renumber face values:
        vertex_mapping = self.renumber_vertices(faces=faces)
        reference_point = np.mean(vertices, axis=0)

        # Begin Writing:
        with open(writeDir, 'w') as file:
            # Writing the header portion of the PLY compliant file:
            file.write('ply\n')
            file.write('format ascii 1.0\n')
            file.write(f'element vertex {vertices.shape[0]}\n')
            file.write('property float32 x\n')
            file.write('property float32 y\n')
            file.write('property float32 z\n')
            file.write(f'element face {faces.shape[0]}\n')
            file.write('property list uint8 int32 vertex_indices\n')
            file.write('end_header\n')

            # Now writing the indices:
            for row in vertices:
                # Round all vertices to 6 decimal places
                file.write(" ".join("{:.6f}".format(x) for x in row) + "\n")

            # And writing the faces:
            for row in faces:
                # We first use the mapping dictionary to ensure we are accessing the correct index positions:
                new_row = []
                for ind in list(row):
                    new_row.append(vertex_mapping[ind])
                # Before writing the face, we need to process it to verify that it is "correct":
                face_vertices = vertices[np.array(new_row)]
                face_normal = self.calculate_face_normal(face_vertices)

                if not self.is_face_clockwise(face_vertices, face_normal, reference_point):
                    # If the face is clockwise, we re-orient properly depending on if rectangular / triangular face:
                    if len(new_row) == 4:
                        new_row = [new_row[0], new_row[3], new_row[2], new_row[1]]
                    else:
                        new_row = [new_row[0], new_row[2], new_row[1]]

                # First we inset into the first index the length of the array for ply compliance. Note that these faces
                # are already "aligned" to be normal-vector facing outwards:
                write_row = np.insert(new_row, 0, len(new_row))
                file.write(" ".join(str(x) for x in write_row) + " \n")

        ## FIXING DATA --> These functions fix the concave issues faces with writing PLY files the way I have.
        my_mesh = trimesh.load_mesh(writeDir)
        my_mesh.fix_normals()

        export = trimesh.exchange.ply.export_ply(mesh=my_mesh, encoding='ascii')
        with open(writeDir, "wb") as f:
            f.write(export)

        # Write filepath2 if needed:
        if filepath is not None and filepath2 is None:
            warnings.warn('You have not specified a filepath 2 value and therefore the script will not write out the '
                          'analysis files for the distance calculator. You must specify a savepath for filepath2 using '
                          'the .aj extension.')
        else:
            fin_verts = []
            for node in self.comparison_graph.nodes():
                pt = self.xyz_from_node(node=node)
                fin_verts.append(pt)
            with open(writeDir2, 'w') as file:
                file.write('Custom binary datafile for analyzing positional data of a generated DNA Nanostructure\n')
                for row in fin_verts:
                    # Round all vertices to 6 decimal places
                    file.write(" ".join("{:.6f}".format(x) for x in row) + "\n")

    def dump_graph_data(self, savepath: str = None) -> None:
        """
        This function is called to dump the design graph data so I can read it in later to debug things.
        :return: Nothing, just exports a file.
        """
        if savepath is None:
            savepath = os.path.join(os.getcwd(), 'SavedGraph')
            if not os.path.exists(savepath):
                os.mkdir(savepath)

        output_file = os.path.join(savepath, 'graph_data_export.pickle')
        copy_of_data = deepcopy(self)
        with open(output_file, 'wb') as f:
            pickle.dump(copy_of_data, f)

    def read_graph_data(self, path_to_pickle: str = None):
        """
        This function simply loads a point cloud so I can validate its point data:
        :return:
        """
        if path_to_pickle is None:
            save_path = os.path.join(os.getcwd(), 'SavedGraph')
            input_file = os.path.join(save_path, 'graph_data_export.pickle')
        else:
            input_file = path_to_pickle
        with open(input_file, 'rb') as f:
            loaded_obj = pickle.load(f)
        return loaded_obj

    def write_out_geometry_json(self, objective_used: str, total_time: float, obj_val: float, savepath: str = None,
                                filename: str = None) -> None:
        """
        This function will write out specific measures for the analysis of the effects of the objective fucntion.
        These will be written to a JSON File so that
        :param objective_used: What objective function was used to store this information for fair comparisons
        :param total_time: Total time of the generative process
        :param obj_val: Final objective function value to record to record against the time
        :param savepath: Where we are saving the output JSON to
        :param filename: What the filename is
        """
        global numDecimals
        if savepath is None:
            savepath = os.path.join(os.getcwd(), 'SavedGraph')
            if not os.path.exists(savepath):
                os.mkdir(savepath)

        if objective_used == 'temp':
            buckling_weight = 8
            C_const = 1 / 8
            regularization = 15
            F_buck, F_rep, N_f, N_e = self.minimize_RBD_forces()
            obj_val = (buckling_weight * F_buck / N_e) + (C_const * F_rep) + (regularization / N_f)

        # Obtain data we need
        mesh_file = self.invalid_mesh(return_mesh=True)

        ## EDGE LENGTHS:
        edge_lengths = []
        for u, v in self.design_graph.edges:
            p1, p2 = self.xyz_from_node(node=u), self.xyz_from_node(node=v)
            distance = self.calculate_3D_distance(p1=p1, p2=p2)
            edge_lengths.append(np.round(distance, numDecimals))
        ## FACE AREAS:
        total_face_areas = list(mesh_file.area_faces)
        numFaces = len(total_face_areas)
        ## FACE ANGLES:
        face_angles = []
        for angles in mesh_file.face_angles:
            for i in angles:
                face_angles.append(np.rad2deg(i))

        ## MESH VOLUMES:
        mesh_volume = round(mesh_file.volume, numDecimals)

        data_to_write = {"Objective Function": objective_used,
                         "Time for Generation (minutes)": round(total_time, numDecimals),
                         "Objective Function Value": round(obj_val, numDecimals),
                         "Edge Lengths": edge_lengths,
                         "Total Number of Faces": numFaces,
                         "Total Face Areas": total_face_areas,
                         "Face Angles": face_angles,
                         "Mesh Volume": mesh_volume}

        # Create the outputfile otherwise:
        if filename is None:
            output_file = os.path.join(savepath, 'geometry_analysis.json')
        else:
            output_file = os.path.join(savepath, filename)
        with open(output_file, 'w') as json_file:
            json.dump(data_to_write, json_file)



def aj_main():
    ## OCTA FRAME:
    """binders = [
        BindingRegion(v1=np.array([25, 25, 0])),
        BindingRegion(v1=np.array([25, 25, 50])),
        BindingRegion(v1=np.array([0, 25, 25])),
        BindingRegion(v1=np.array([50, 25, 25])),
        BindingRegion(v1=np.array([25, 0, 25])),
        BindingRegion(v1=np.array([25, 50, 25]))
    ]"""

    ### Asymmetric:
    binders = [
        BindingRegion(v1=np.array([25, 0, 25])),
        BindingRegion(v1=np.array([15, 50, 25]), v2=np.array([35, 50, 25])),
        BindingRegion(v1=np.array([35, 35, 50]), v2=np.array([35, 15, 50]),
                      v3=np.array([15, 15, 50]), v4=np.array([15, 35, 50])),
        BindingRegion(v1=np.array([25, 25, 0]))
    ]

    ### Beam:
    """cross_section_edge_length = 37
    beam_length = 70
    maxY = cross_section_edge_length + 15  # Give a "large" space for exploration!
    maxZ = cross_section_edge_length + 15
    sq_by_2 = cross_section_edge_length / 2
    binders = [
        BindingRegion(v1=np.array([0, (maxY / 2) - sq_by_2, (maxZ / 2) - sq_by_2]),
                          v2=np.array([0, (maxY / 2) - sq_by_2, (maxZ / 2) + sq_by_2]),
                          v3=np.array([0, (maxY / 2) + sq_by_2, (maxZ / 2) + sq_by_2]),
                          v4=np.array([0, (maxY / 2) + sq_by_2, (maxZ / 2) - sq_by_2])),
        BindingRegion(v1=np.array([beam_length, (maxY / 2) - sq_by_2, (maxZ / 2) - sq_by_2]),
                          v2=np.array([beam_length, (maxY / 2) - sq_by_2, (maxZ / 2) + sq_by_2]),
                          v3=np.array([beam_length, (maxY / 2) + sq_by_2, (maxZ / 2) + sq_by_2]),
                          v4=np.array([beam_length, (maxY / 2) + sq_by_2, (maxZ / 2) - sq_by_2]))
    ]"""

    ## NP
    NPs = [
        NPBox(c1=np.array([0, 15, 15]), c2=np.array([20, 35, 35]))
    ]
    X = 50
    Y = 50
    Z = 50
    design_space = DesignSpace3D(maxX=X, maxY=Y, maxZ=Z, NP_Regions=NPs, binding_regions=binders,
                                 min_edge_length=31, routing_software='DAEDALUS', debug_mode=False)
    # FOR BEAM:
    """design_space = DesignSpace3D(maxX=beam_length, maxY=maxY, maxZ=maxZ, NP_Regions=[], binding_regions=binders,
                                 min_edge_length=31, routing_software='DAEDALUS', debug_mode=False)"""


    design_space.generate_start_shape()
    #design_space.write_out_geometry_json(objective_used='temp', total_time=5.5, obj_val = 1)


    """design_space.generate_start_shape()
    design_space.display_graph()
    design_space.retriangulate_face_rule()
    design_space.move_vert_in_average_normal_dir_rule()
    design_space.dump_graph_data()
    del design``_space
    design_space = DesignSpace3D(maxX=X, maxY=Y, maxZ=Z, NP_Regions=[], binding_regions=binders,
                                 min_edge_length=31, routing_software='vHelix', debug_mode=True)
    #design_space.read_graph_data()
    design_space = design_space.read_graph_data()
    design_space.display_graph()"""

    design_space.display_graph(use_in_dash=True)
    design_space.display_graph_as_cylindrical_rep()
    #design_space.retriangulate_face_rule()
    #temp = design_space.calculate_cost_function()




if __name__ == '__main__':
    aj_main()
