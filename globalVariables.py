import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import copy
import math
from constants import *

class grid_class:
    def __init__(self):
        self.state = np.random.randint(2, states, size = (n,n))
        self.dislocation_densities = [[P_initial for l in range(n)] for m in range(n)]
        self.dynamic_recrystallization_number = np.zeros((n, n))
        self.orientation = [[np.random.randint(1, 360, size = (n,n)), np.random.randint(1, 180, size = (n,n)), np.random.randint(1, 360, size = (n,n))]]


current_strain = 0.0172  
true_strain = []
true_stress = []

global new_grains
new_grains = 0

global total_grains
total_grains = n * n

global state_to_dislocation
state_to_dislocation = []

global misorientations
misorientations = []

global delta_dislocations
delta_dislocations = []

global velocities
velocities = []

global state_to_orientation_ph
state_to_orientation_ph = []

global state_to_orientation_th
state_to_orientation_th = []

global state_to_orientation_om
state_to_orientation_om = []


def interpolate_array(array):
    current_length = len(array)
    for i in range(current_length - 1):
        array.insert(2 * i + 1, ((array[2 * i] + array[2 * i + 1]) / 2))
    return array

strains = [0.017232,
0.065399,
0.14911,
0.234872,
0.301806,
0.38547,
0.450306,
0.538158,
0.600907,
0.686673,
0.751526,
0.833113,
0.902147,
0.975369,
1.0465,
1.11972,
1.19712]

strains = interpolate_array(strains)

for g in range(10):
    strains = interpolate_array(strains)

iterations = len(strains)
global k
k = 0

border_grid = np.zeros((n, n))
recrystallized_grid = np.zeros((n, n))
cell_strain = np.zeros((n, n))
original_grid = grid_class()
updation_grid = copy.deepcopy(original_grid)