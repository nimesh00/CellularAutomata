#! /usr/bin/env python3

import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import copy
import math

states = 50
n = 100
iterations = 30
strain_rate = 0.01
P_initial = 6.022 * (10 ** 14)
P_RX_initial = 6 * (10 ** 13)
P_max = 8.08 * (10 ** 14)
P_cr = 8.08 * (10 ** 14)
e_cr = 0.15
critical_misorientation = 15
k1 = 7.78 * 10 ** 8
k2 = 27.09
cd = 2 * 10 ** -6
Mo = 1.4 * (10 ** -5)
Tm = 1453  
mu_o = 7.89 * 10 ** 4
Current_Temp = 1313
b = 2.49 * 10 ** -10
alpha = 0.5
Poisson_ratio = 0.3

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


def grain_initialization():
    iterations_init = 100000
    mat = np.zeros((n,n), dtype=int)
    for m in range(states):
        i = int(random.randrange(n))
        j = int(random.randrange(n))
        mat[i][j] = 1 + random.randrange(states)

    for m in range(iterations_init):
        i = random.randint(0,n) - 1
        j = random.randint(0,n) - 1
        if m % 2 == 0:
            neighbour_indices = [[i - 1, j - 1], [i - 1, j], 
                                [i, j - 1],                      [i, j + 1],
                                                [i + 1, j], [i + 1, j + 1]]
        else:
            neighbour_indices = [                [i - 1, j], [i - 1, j + 1],
                                [i, j - 1],                      [i, j + 1],
                                [i + 1, j - 1], [i + 1, j],               ]

    
        if mat[i][j] != 0:
            for indices in neighbour_indices:
                if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
                    continue

                elif mat[indices[0]][indices[1]] == 0:
                    mat[indices[0]][indices[1]] = mat[i][j]

                else:
                    continue
    
    return mat


def most_common(List): 
    counter = 0
    num = List[0] 
      
    for i in List: 
        curr_frequency = List.count(i) 
        if(curr_frequency> counter): 
            counter = curr_frequency 
            num = i 
  
    return num


def monte_carlo_initialization(mat):
    states_list = []
    iterations_mc = 100000
    for i in range(states):
        states_list.append(i+1)    

    for m in range(iterations_mc):
        i = random.randint(1,n) - 1
        j = random.randint(1,n) - 1
        if m % 2 == 0:
            neighbour_indices = [[i - 1, j - 1], [i - 1, j], 
                                [i, j - 1],                      [i, j + 1],
                                                [i + 1, j], [i + 1, j + 1]]
        else:
            neighbour_indices = [                [i - 1, j], [i - 1, j + 1],
                                [i, j - 1],                      [i, j + 1],
                                [i + 1, j - 1], [i + 1, j],               ]

        not_similar_i = 0
        not_similar_f = 0
        switch = 0
        original = mat[i][j]
        neighbourhood = []

        #find the initial energy before switching 
        for indices in neighbour_indices:
            if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
                continue
            neighbourhood.append(mat[indices[0]][indices[1]])
            if mat[indices[0]][indices[1]] == mat[i][j]:
                continue
            else: 
                not_similar_i += 1
                #this would give me the total number of different cells in neighbour initially (or initial energy Ei)

        switch = most_common(neighbourhood)
        #gives the element I would switch and calculate the energy with

        if switch == 0:
            mat[i][j] = original
        else:
            mat[i][j] = switch
    
        #find the final energy after switching         
        for indices in neighbour_indices:
            if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
                continue
            neighbourhood.append(mat[indices[0]][indices[1]])
            if mat[indices[0]][indices[1]] == mat[i][j]:
                continue
            else: 
                not_similar_f += 1   #this would give me the total number of different cells in neighbour finally (or final energy Ef)


        #calculating whether switch must be accepted or not.

        delta_e = not_similar_f - not_similar_i

        if delta_e <= 0:
            mat[i][j] = switch
        else:
            prob = abs(delta_e)/8
            rand = random.random()
            if prob >= rand:
                mat[i][j] = switch
            else:
                mat[i][j] = original

    return mat

def mu(temp):
    mu = mu_o * (1 - 0.64 * (temp - 27) / (Tm + 273))
    return mu


def mis_orientation(ph_1, th_1, om_1, ph_2, th_2, om_2):
    #this function calculates and returns misorientation between two adjacent cells
    #here ph is 'phi', th is 'theta', om is 'omega'
    R_1 = [[math.cos(om_1) * math.cos(ph_1) - math.sin(ph_1) * math.sin(om_1) * math.cos(th_1),    math.cos(om_1) * math.sin(ph_1) + math.cos(ph_1) * math.sin(om_1) * math.cos(th_1),   math.sin(th_1) * math.sin(om_1)],
           [-math.sin(om_1) * math.cos(ph_1) - math.sin(ph_1) * math.cos(om_1) * math.cos(th_1),  -math.sin(ph_1) * math.sin(om_1) + math.cos(ph_1) * math.cos(om_1) * math.cos(th_1),   math.cos(om_1) * math.sin(th_1)],
           [math.sin(ph_1) * math.sin(th_1),                                                      -math.sin(th_1) * math.cos(ph_1),                                                                       math.cos(th_1)]]

    R_2 = [[math.cos(om_2) * math.cos(ph_2) - math.sin(ph_2) * math.sin(om_2) * math.cos(th_2),    math.cos(om_2) * math.sin(ph_2) + math.cos(ph_2) * math.sin(om_2) * math.cos(th_2),   math.sin(th_2) * math.sin(om_2)],
           [-math.sin(om_2) * math.cos(ph_2) - math.sin(ph_2) * math.cos(om_2) * math.cos(th_2),  -math.sin(ph_2) * math.sin(om_2) + math.cos(ph_2) * math.cos(om_2) * math.cos(th_2),   math.cos(om_2) * math.sin(th_2)],
           [math.sin(ph_2) * math.sin(th_2),                                                      -math.sin(th_2) * math.cos(ph_2),                                                                       math.cos(th_2)]]

    trans_R_2 = [[R_2[j][i] for j in range(len(R_2))] for i in range(len(R_2[0]))]

    multiply = [[0,0,0],  
                [0,0,0],  
                [0,0,0]]  
   
    for i in range(len(R_1)):  
       for j in range(len(trans_R_2[0])):  
           for k in range(len(trans_R_2)):  
               multiply[i][j] += R_1[i][k] * trans_R_2[k][j]  
    
    trace = 0
    for i in range(3):
        for j in range(3):
            if i == j:
                trace += multiply[i][j]

    final = (trace - 1) / 2
    mis_orient = np.arccos(final) * 180 / math.pi
    return mis_orient

    
def grain_velocity(delta_dislocation_btw_adj_cells, misorientation_btw_adj_cells):                          
    mobility  = Mo * (1 - np.exp( -5 * ((float(misorientation_btw_adj_cells+1) / critical_misorientation) ** 4)))

    gama_o = (mu(Current_Temp) * b * critical_misorientation) / (4 * math.pi * (1 - Poisson_ratio))
    if misorientation_btw_adj_cells >= critical_misorientation:
        gama = gama_o                                                   #gama is grain_boundary_energy 
    elif misorientation_btw_adj_cells == 0:
        gama = 0
    else:
        gama = gama_o * (misorientation_btw_adj_cells / critical_misorientation) * (1 - np.log(misorientation_btw_adj_cells / critical_misorientation))

    tau = 0.5 * mu(Current_Temp) * (b ** 2)                             #tau is dislocation_line_energy
    r = 1 #np.sqrt(200 * cd * cd / math.pi)
    P = (tau * delta_dislocation_btw_adj_cells) - (2 * gama / r)                      #P is stored_deformation_energy

    grain_velo = mobility * P
    return grain_velo


def strain(avg_velocity):                                                                                                                                                                                                                                        
    delta_t = cd / avg_velocity
    global current_strain
    delta_strain = strain_rate * delta_t
    current_strain +=  delta_strain
    return current_strain, delta_strain


def update_dislocation_density(P_prev_time_step, delta_strn):                                                                                                                                    
    delta_p_for_same_cell = (k1 * np.sqrt(P_prev_time_step) - k2 * P_prev_time_step) * (delta_strn)
    P_new = P_prev_time_step + delta_p_for_same_cell
    return P_new


def update_recrystallized_cell_state(i, j):
    nucleation_probability = original_grid.dislocation_densities[i][j] / np.max(original_grid.dislocation_densities)

    random_number = np.random.random()
    
    if random_number <= nucleation_probability:
        updation_grid.state[i][j] = np.random.randint(2, states)
        updation_grid.dislocation_densities[i][j] = P_RX_initial
        updation_grid.dynamic_recrystallization_number[i][j] = 1
        updation_grid.orientation[0][0][i][j] = state_to_orientation_ph[original_grid.state[i][j]]                                               
        updation_grid.orientation[0][1][i][j] = state_to_orientation_th[original_grid.state[i][j]]
        updation_grid.orientation[0][2][i][j] = state_to_orientation_om[original_grid.state[i][j]]
        #print(i, j, " nucleated")
        return 1
    return 0


def simulated_stress():
    global total_grains
    average_dislocation_density = np.mean(original_grid.dislocation_densities)
    simulated_stress = alpha * mu(Current_Temp) * b * np.sqrt(average_dislocation_density)
    return simulated_stress


def cell_on_border(i, j):
    #neighbour_indices = [                [i - 1, j],               
    #                     [i, j - 1],                     [i, j + 1],
    #                                     [i + 1, j],               ]

    neighbour_indices = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],               
                         [i, j - 1],                     [i, j + 1],
                         [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]
    
    for indices in neighbour_indices:
        if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
            continue
        if original_grid.state[i][j] != original_grid.state[indices[0]][indices[1]]:
            border_grid[i][j] = 1
            return True
    return False


def check_mat_for_zero(matrix):
    height, width = matrix.shape

    for i in range(height):
        for j in range(width):
            if matrix[i][j] == 0:
                return True
    return False


def near_recrystallized_cell(i, j):
    #neighbour_indices = [                [i - 1, j],               
    #                     [i, j - 1],                     [i, j + 1],
    #                                     [i + 1, j],               ]
    
    neighbour_indices = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],               
                         [i, j - 1],                     [i, j + 1],
                         [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]
    
    for indices in neighbour_indices:
        if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
            continue
        if original_grid.dynamic_recrystallization_number[indices[0]][indices[1]] == 1:
            recrystallized_grid[i][j] = 1
            return True
    return False


def propagateGrainBoundary(i, j, k):
    if k % 2 == 0:
        neighbour_indices = [[i - 1, j - 1], [i - 1, j], 
                            [i, j - 1],                      [i, j + 1],
                                            [i + 1, j], [i + 1, j + 1]]
    else:
        neighbour_indices = [                [i - 1, j], [i - 1, j + 1],
                            [i, j - 1],                      [i, j + 1],
                            [i + 1, j - 1], [i + 1, j],               ]
    
    favoured_indices = [0, 0]
    found_nucleus = False
    min_p_neighbour = np.max(original_grid.dislocation_densities)
    for indices in neighbour_indices:
        if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
            continue
        
        if original_grid.dynamic_recrystallization_number[indices[0]][indices[1]] == 1:
            found_nucleus = True
            p_neighbour = original_grid.dislocation_densities[indices[0]][indices[1]]
            if p_neighbour < min_p_neighbour:
                min_p_neighbour = p_neighbour
                favoured_indices = indices

    if found_nucleus is False:
        return 0

    propagation_probability = original_grid.dislocation_densities[i][j] * ((float(k2) / k1) ** 2)
    random_number = np.random.random()

    if random_number <= propagation_probability:
        updation_grid.orientation[0][0][i][j] =                original_grid.orientation[0][0][favoured_indices[0]][favoured_indices[1]]
        updation_grid.orientation[0][1][i][j] =                original_grid.orientation[0][1][favoured_indices[0]][favoured_indices[1]]
        updation_grid.orientation[0][2][i][j] =                original_grid.orientation[0][2][favoured_indices[0]][favoured_indices[1]]
        updation_grid.state[i][j] =                            original_grid.state[favoured_indices[0]][favoured_indices[1]]
        updation_grid.dislocation_densities[i][j] =            P_RX_initial
        updation_grid.dynamic_recrystallization_number[i][j] = 1
        return 1
    return 0


def update_grid_state(plt, fig , ax, cmap, bounds):
    N = states
    x = np.zeros((n*n))
    y = np.zeros((n*n))
    tag = np.zeros(n*n)
    
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    for i in range(n):
        for j in range(n):
            x[i*n + j] = i
            y[i*n + j] = j
            tag[i*n + j] = original_grid.state[i][j]
    
    scat = ax.scatter(x, y, marker='s', c = tag, cmap = cmap, norm=norm)

    return scat

def initialize_dislocation_density():
    possible_random_numbers = np.random.random((states+1))
    dislocation_densities = [[P_initial for i in range(n)] for j in range(n)]                                       #creates states-2 diff random values only?
    for i in range(n):
        for j in range(n):
            dislocation_densities[i][j] = (2 * possible_random_numbers[original_grid.state[i][j]] + 6) * (10 ** 14)
            #dislocation_densities[i][j] = (2 * np.random.random() + 6) * (10 ** 14)
    return dislocation_densities


def initialize_orientation():
    global state_to_orientation_ph
    global state_to_orientation_th
    global state_to_orientation_om
    
    state_to_orientation_ph = np.random.randint(1, 360, size = (states + 1))
    state_to_orientation_th = np.random.randint(1, 180, size = (states + 1))
    state_to_orientation_om = np.random.randint(1, 360, size = (states + 1))
    original_grid.orientation[0][0] = [[0 for i in range(n)] for j in range(n)]
    original_grid.orientation[0][1] = [[0 for i in range(n)] for j in range(n)]
    original_grid.orientation[0][2] = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            original_grid.orientation[0][0][i][j] = state_to_orientation_ph[original_grid.state[i][j]]
            original_grid.orientation[0][1][i][j] = state_to_orientation_th[original_grid.state[i][j]]
            original_grid.orientation[0][2][i][j] = state_to_orientation_om[original_grid.state[i][j]]
    return original_grid.orientation

def potential_nucleus(i, j):
    if original_grid.dynamic_recrystallization_number[i][j] == 1:
        return False
    #neighbour_indices = [                [i - 1, j],               
    #                     [i, j - 1],                     [i, j + 1],
    #                                     [i + 1, j],               ]
    
    neighbour_indices = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],               
                         [i, j - 1],                     [i, j + 1],
                         [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]
    
    #critical_dP = 0.001 * np.max(original_grid.dislocation_densities)

    for indices in neighbour_indices:
        if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
            continue
        elif (original_grid.orientation[0][0][i][j] == original_grid.orientation[0][0][indices[0]][indices[1]]) and \
             (original_grid.orientation[0][1][i][j] == original_grid.orientation[0][1][indices[0]][indices[1]]) and \
             (original_grid.orientation[0][2][i][j] == original_grid.orientation[0][2][indices[0]][indices[1]]):
            continue
        misorientation = mis_orientation(original_grid.orientation[0][0][i][j], 
                                        original_grid.orientation[0][1][i][j], 
                                        original_grid.orientation[0][2][i][j], 
                                        original_grid.orientation[0][0][indices[0]][indices[1]], 
                                        original_grid.orientation[0][1][indices[0]][indices[1]], 
                                        original_grid.orientation[0][2][indices[0]][indices[1]])

        condition1 = (misorientation >= critical_misorientation)

        #dP = original_grid.dislocation_densities[i][j] - original_grid.dislocation_densities[indices[0]][indices[1]]
        condition2 = (original_grid.dislocation_densities[i][j] >= P_cr)

        # print("misorientation, dP: ", misorientation, (dP / critical_dP))
        # print("cell, neighbour: ", original_grid.dislocation_densities[i][j], original_grid.dislocation_densities[indices[0]][indices[1]])
        # print("min, max: ", np.min(original_grid.dislocation_densities), np.max(original_grid.dislocation_densities))

        # print(condition1, condition2)
        
        if condition1 and condition2:
            return True
    return False

def potential_growth(i, j, k):
    if original_grid.dynamic_recrystallization_number[i][j] == 1:
        return False
    if k % 2 == 0:
        neighbour_indices = [[i - 1, j - 1], [i - 1, j], 
                            [i, j - 1],                      [i, j + 1],
                                            [i + 1, j], [i + 1, j + 1]]
    else:
        neighbour_indices = [                [i - 1, j], [i - 1, j + 1],
                            [i, j - 1],                      [i, j + 1],
                            [i + 1, j - 1], [i + 1, j],               ]

    #critical_dP = 0.001 * np.max(original_grid.dislocation_densities)

    for indices in neighbour_indices:
        if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
            continue
        elif (original_grid.orientation[0][0][i][j] == original_grid.orientation[0][0][indices[0]][indices[1]]) and \
             (original_grid.orientation[0][1][i][j] == original_grid.orientation[0][1][indices[0]][indices[1]]) and \
             (original_grid.orientation[0][2][i][j] == original_grid.orientation[0][2][indices[0]][indices[1]]):
            continue
        misorientation = mis_orientation(original_grid.orientation[0][0][i][j], 
                                        original_grid.orientation[0][1][i][j], 
                                        original_grid.orientation[0][2][i][j], 
                                        original_grid.orientation[0][0][indices[0]][indices[1]], 
                                        original_grid.orientation[0][1][indices[0]][indices[1]], 
                                        original_grid.orientation[0][2][indices[0]][indices[1]])

        condition1 = (misorientation >= critical_misorientation)

        if original_grid.dynamic_recrystallization_number[indices[0]][indices[1]] == 1 and condition1:
            return True
    return False


def main():
    global true_strain
    global true_stress
    global total_grains
    global original_grid
    global updation_grid
    global misorientations
    global delta_dislocations
    global velocities

    original_grid.state = monte_carlo_initialization(grain_initialization())
    original_grid.dislocation_densities = initialize_dislocation_density()
    original_grid.orientation = initialize_orientation()
    
    k = 0
    N = states
    xDrx = []
    
    # setup the plot
    fig, ax = plt.subplots(1,3, figsize=(10,10))
    cmap = plt.cm.jet
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    cmap = plt.cm.jet
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    bounds = np.linspace(0,N,N+1)

    cb = plt.colorbar(update_grid_state(plt, fig, ax[0], cmap, bounds), spacing='proportional',ticks=bounds)
    cb.set_label('Custom cbar')
    

    update_grid_state(plt, fig, ax[0], cmap, bounds)
    plt.pause(0.01)

    N_c = n * n
    N_r = 0
    N_unchanged_grains = n * n
    ret = 0
    nucleatinon_step = True
    nucleation_happened = False

    update_indices = []
    propagate_indices = []

    updation_grid = copy.deepcopy(original_grid)
    
    start_process = False

    while (k < iterations):
        print("min, max: ", np.min(original_grid.dislocation_densities), np.max(original_grid.dislocation_densities))
        if np.sum(original_grid.dynamic_recrystallization_number) == n * n:
            print("Bas aaj k liye itna hi!!")
            break
        #     print("All cells recrystallized!! Using the new grid and starting the process again!!")
        #     updation_grid.dynamic_recrystallization_number = np.zeros((n, n))
        prev_grid = copy.deepcopy(original_grid)
        counter_propagate1 = 0
        counter_propagate2 = 0
        #for the overall strain calculation
        for i in range(n):
            for j in range(n):
                neighbour_indices = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],
                                    [i, j - 1],                      [i, j + 1],
                                    [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]
            
                if cell_on_border(i, j):
                    for indices in neighbour_indices:
                        if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
                            continue
                        #elif original_grid.state[i][j] == original_grid.state[indices[0]][indices[1]]:
                        elif (original_grid.orientation[0][0][i][j] == original_grid.orientation[0][0][indices[0]][indices[1]]) and (original_grid.orientation[0][1][i][j] == original_grid.orientation[0][1][indices[0]][indices[1]]) and (original_grid.orientation[0][2][i][j] == original_grid.orientation[0][2][indices[0]][indices[1]]):
                            continue
                        else:
                            misorientations += [mis_orientation(original_grid.orientation[0][0][i][j], 
                                                                original_grid.orientation[0][1][i][j], 
                                                                original_grid.orientation[0][2][i][j], 
                                                                original_grid.orientation[0][0][indices[0]][indices[1]], 
                                                                original_grid.orientation[0][1][indices[0]][indices[1]], 
                                                                original_grid.orientation[0][2][indices[0]][indices[1]])]
                            delta_dislocations += [abs(original_grid.dislocation_densities[i][j] - original_grid.dislocation_densities[indices[0]][indices[1]])]

        for i in range(len(misorientations)):
            velocities += [grain_velocity(delta_dislocations[i], misorientations[i])]
            
        stn, delta_strain = strain(np.mean(velocities))

        avg_dis_den = np.mean(original_grid.dislocation_densities)
        if stn >= e_cr:
            start_process = True
        #your normal code continues from here
        for i in range(n):
            for j in range(n):
                if original_grid.dynamic_recrystallization_number[i][j] == 1:
                    # updation_grid.dislocation_densities[i][j] = update_dislocation_density(updation_grid.dislocation_densities[i][j], delta_strain)
                    continue
                ret = 0
                cell_on_border(i, j)
                near_recrystallized_cell(i, j)
    
                avg_dis_den = np.mean(original_grid.dislocation_densities)
                if stn >= e_cr:
                    start_process = True
                # if avg_dis_den >= P_cr:
                if start_process is True:
                    # print("Upar dekho")
                    if potential_growth(i, j, k):
                        ret = propagateGrainBoundary(i, j, k)
                    if ret == 1:
                        continue
                    if potential_nucleus(i, j):
                        ret = update_recrystallized_cell_state(i, j)
                    if ret == 1:
                        # print("Nucleate: ", i, j)
                        continue
                    # print("Niche dekho")
                    updation_grid.dislocation_densities[i][j] = update_dislocation_density(updation_grid.dislocation_densities[i][j], delta_strain)
                else:
                    updation_grid.dislocation_densities[i][j] = update_dislocation_density(updation_grid.dislocation_densities[i][j], delta_strain)

        if nucleation_happened:
            nucleatinon_step = False
        original_grid = copy.deepcopy(updation_grid)

        unq, cnts = np.unique(original_grid.dynamic_recrystallization_number, return_counts=True)
        unq_array = dict(zip(unq, cnts))
        #print("updated cells: ", unq_array)
        difference_matrix = original_grid.dynamic_recrystallization_number - prev_grid.dynamic_recrystallization_number
        # N_unchanged_grains = difference_matrix.count(0)
        unique, counts = np.unique(difference_matrix, return_counts=True)
        unique_array = dict(zip(unique, counts))
        # print(unique_array)
        N_unchanged_grains = unique_array[0]

        # total number of grains
        global new_grains
        new_grains = N_c - N_unchanged_grains
        N_r += N_c - N_unchanged_grains                                                                                    
        
        total_grains = N_c

        # print("Total Grains: ", total_grains)
        
        xDrx = xDrx + [float(N_r) / N_c]
        # print(xDrx)

        net_strain, _ = strain(np.mean(velocities))
        #net_strain = np.mean(cell_strain)
        net_stress = simulated_stress()
        #print("Net strain value:", net_strain)
        # print("Net stress value:", net_stress)
        # true_strain = true_strain + [net_strain]
        true_strain += [net_strain] 
        true_stress += [net_stress]
        #print("True Strain: ", true_strain)
        print("Iteration",k, "completed")

        # print(grid)
        # print(dislocation_densities)
        # pw.plot(true_strain, true_stress, pen='r')
        k += 1
        # print(grid)
        # BLOCK 1 START
        update_grid_state(plt, fig, ax[0], cmap, bounds)
        ax[1].plot(true_strain, true_stress, 'b-')
        ax[2].plot(true_strain, xDrx, 'r-')
        plt.pause(0.01)        
        # BLOCK 1 END
    # BLOCK 2 START
    # update_grid_state(plt, fig, ax[0], cmap, bounds)
    # ax[1].plot(true_strain, true_stress, 'b-')
    # ax[2].plot(true_strain, xDrx, 'r-')
    #BLOCK 2 END
    plt.show()

if __name__ == "__main__":
    main()















