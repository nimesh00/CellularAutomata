'''
    COMMENT BLOCK 1 and UNCOMMENT BLOCK 2 : Updations first and then end plot. End results only
    UNCOMMENT BLOCK 1 and COMMENT BLOCK 2: Live updations after each iterations. Slower(can be optimized)
    Enter the correct values of the constants
'''

import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import copy
import math
# from helper import grid_class

states = 20
n = 100
iterations = 30
strain_rate = 0.01
P_initial = 6.022 * (10 ** 14)
P_max = 8.214 * (10 ** 14)
P_cr = 8.214 * (10 ** 14)
critical_orientation = 15
k1 = 7.78 * 10 ** 8
k2 = 27.09
cd = 2 * 10 ** -9
# Mo = 1 * (10 ** 5)
Mo = 1.4 * (10 ** -11)
Tm = 1453  # degree kelvin
mu_o = 7.89 * 10 ** 4
Current_Temp = 1313
b = 2.49 * 10 ** -10
alpha = 0.5

class grid_class:
    def __init__(self):
        self.state = np.random.randint(2, states, size = (n,n))
        self.dislocation_densities = [[P_initial for l in range(n)] for m in range(n)]
        self.dynamic_recrystallization_number = np.zeros((n, n))
        self.orientation = np.random.randint(1, 180, size = (n,n))

current_strain = 0
true_strain = []
true_stress = []
global new_grains
new_grains = 0

global total_grains
total_grains = n * n

global state_to_orientation
state_to_orientation = []

def interpolate_array(array):
    #doubles the length
    current_length = len(array)
    # print("original array: ", array)
    for i in range(current_length - 1):
        # print(i, " elements added")
        # print("New length: ", len(array))
        array.insert(2 * i + 1, ((array[2 * i] + array[2 * i + 1]) / 2))
        # print("new array: ", array)
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

for g in range(11):
    strains = interpolate_array(strains)
    print("New length: ", len(strains))

iterations = len(strains)
global k
k = 0

border_grid = np.zeros((n, n))
recrystallized_grid = np.zeros((n, n))
# orientation = np.random.randint(1, 180, size = (n,n))
# grid = np.random.randint(2, states, size = (n,n))
# updated_grid = copy.copy(grid)
# dislocation_densities = [[P_initial for l in range(n)] for m in range(n)]
# dynamic_recrystallization_number = np.zeros((n, n))
cell_strain = np.zeros((n, n))

original_grid = grid_class()
updation_grid = copy.deepcopy(original_grid)


def mu(temp):
    mu = mu_o * (1 - 0.64 * (temp - 27) / (Tm + 273))
    return mu

def simulated_stress():
    # average_dislocation_energy = np.mean(dislocation_densities)
    global total_grains
    average_dislocation_energy = np.mean(original_grid.dislocation_densities)
    # print("Densities: ", average_dislocation_energy, P_cr)
    simulated_stress = alpha * mu(Current_Temp) * b * np.sqrt(average_dislocation_energy)
    return simulated_stress

def mobility(orientation):
    return (Mo * (1 - np.exp( -5 * ((float(orientation) / critical_orientation) ** 4))))

def strain(orientation):
    # print('Cd: ', cd)
    # print('k1: ', k1)
    # print('k2: ', k2)
    # print('Mo: ', Mo)
    # print('orientation: ', orientation)
    # print('critical orientation: ', critical_orientation)
    delta_t = (cd * ((float(k2) / k1) ** 2)) / (0.5 * mu(Current_Temp) * (b ** 2)) * mobility(orientation)
    global current_strain
    current_strain +=  strain_rate * delta_t
    # print("current_strain: ", current_strain)
    return current_strain


def dislocation_energy(P_prev, orientation_current):
    # print('P_previous: ', P_prev)
    global k
    if (k == 0):
        delta_p = (k1 * np.sqrt(P_prev) - k2 * P_prev) * (strains[k])
    else:    
        delta_p = (k1 * np.sqrt(P_prev) - k2 * P_prev) * (strains[k] - strains[k - 1])
    # print('change in energy: ', delta_p)
    P_new = P_prev + delta_p
    return P_new

def update_state(P_prev, orientation_current, state_current):
    P_new = dislocation_energy(P_prev, orientation_current)
    new_state = state_current
    # print("Energy Difference", P_new - P_cr)
    P_max = np.max(original_grid.dislocation_densities)
    choice = 0
    # print("P_new, P_max: ", P_new, P_max)
    nucleation_probability = P_new / P_max
    random_number = np.random.randint(0, 100)
    random_number = float(random_number) / 100
    # print("Random Number, Nucleation Prob", random_number, nucleation_probability)
    if random_number < nucleation_probability:
        choice = 1
        new_state = np.random.randint(2, states)

    return choice, new_state

def update_cell_state(i, j):
    nucleation_probability = original_grid.dislocation_densities[i][j] / np.max(original_grid.dislocation_densities)

    random_number = np.random.random()
    
    if random_number <= nucleation_probability:
        updation_grid.state[i][j] = np.random.randint(2, states)
        updation_grid.dislocation_densities[i][j] = 0
        updation_grid.dynamic_recrystallization_number[i][j] = 1
        updation_grid.orientation[i][j] = state_to_orientation[original_grid.state[i][j]]
        # print(i, j, " nucleated")
        return 1
    return 0

def propagateGrainBoundary(i, j, k):
    # moore's neighbourhood
    # neighbour_indices = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],
    #                     [i, j - 1],                      [i, j + 1],
    #                     [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]
    
    # von-neuman neighbourhood
    #neighbour_indices = [[i - 1, j], [i, j - 1], [i, j + 1], [i + 1, j]]
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
            # delta_p = abs(dislocation_densities[i][j] - dislocation_densities[indices[0]][indices[1]])
            p_neighbour = original_grid.dislocation_densities[indices[0]][indices[1]]
            if p_neighbour < min_p_neighbour:
                min_p_neighbour = p_neighbour
                favoured_indices = indices

    if found_nucleus is False:
        return 0

    # print(i, j, " choosing: ", favoured_indices)

    stored_deformation_energy = 0.5 * original_grid.dislocation_densities[i][j] * mu(Current_Temp) * (b ** 2)
    delta_orientation = abs(original_grid.orientation[i][j] - original_grid.orientation[favoured_indices[0]][favoured_indices[1]])
    # current_mobility = mobility(delta_orientation)
    # vel = current_mobility * stored_deformation_energy
    # vel_max = current_mobility * (0.5 * mu(Current_Temp) * (b ** 2)) * ((k1 / k2) ** 2)
    # nucleation_probability = vel / vel_max
    nucleation_probability = original_grid.dislocation_densities[i][j] * ((float(k2) / k1) ** 2)
    # print("nucleation P: ", nucleation_probability)
    random_number = np.random.random()

    if random_number <= nucleation_probability:
        updation_grid.orientation[i][j] = original_grid.orientation[favoured_indices[0]][favoured_indices[1]]
        updation_grid.state[i][j] = original_grid.state[favoured_indices[0]][favoured_indices[1]]
        updation_grid.dislocation_densities[i][j] = original_grid.dislocation_densities[favoured_indices[0]][favoured_indices[1]]
        updation_grid.dynamic_recrystallization_number[i][j] = 1
        return 1
    return 0

def update_grid_state(plt, fig , ax, cmap, bounds):
    N = states
    x = np.zeros((n*n))
    y = np.zeros((n*n))
    tag = np.zeros(n*n)
    
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    # define the bins and normalize

    for i in range(n):
        for j in range(n):
            x[i*n + j] = i
            y[i*n + j] = j
            tag[i*n + j] = original_grid.state[i][j]
    
    scat = ax.scatter(x, y, marker='s', c = tag, cmap = cmap, norm=norm)

    return scat


# def updating_plot(pw, x, y):
#     pw = pg.plot(x, y, pen = 'r')

def check_mat_for_zero(matrix):
    height, width = matrix.shape

    for i in range(height):
        for j in range(width):
            if matrix[i][j] == 0:
                return True
    return False

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
        #print(m)
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
            
def grain_initialization():
    iterations_init = 100000
    mat = np.zeros((n,n), dtype=int)
    # print(mat)
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

def vornoi_initialization(mat):
    matx, maty = n, n
    nx = []
    ny = []
    ns = []
    for i in range(states):
        nx.append(random.randrange(matx))
        ny.append(random.randrange(maty))
        ns.append(random.randrange(states))
    for y in range(maty):
        for x in range(matx):
            dmin = math.hypot(matx-1, maty-1)
            j = states - 1
            for i in range(states):
                d = math.hypot(nx[i]-x, ny[i]-y)
                if d < dmin:
                    dmin = d
                    j = i
            mat[y][x] = j
    print(mat)
    return mat

def initialize_dislocation_density():
    possible_random_numbers = np.random.random((states))
    dislocation_densities = [[P_initial for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            # dislocation_densities[i][j] = (2 * possible_random_numbers[original_grid.state[i][j]] + 7) * (10 ** 14)
            dislocation_densities[i][j] = (1 * np.random.random() + 6) * (10 ** 14)
    return dislocation_densities

def initialize_orientation():
    global state_to_orientation
    state_to_orientation = np.random.randint(1, 180, size = (states + 1))
    orientation = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            orientation[i][j] = state_to_orientation[original_grid.state[i][j]]
    return orientation

def near_recrystallized_cell(i, j):
    #m moore's neighbourhood
    neighbour_indices = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],
                         [i, j - 1],                     [i, j + 1],
                         [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]
    
    # neighbours = dynamic_recrystallization_number[i-1 : i+2, j-1 : j+2]
    
    # von-neuman neighbourhood
    #neighbour_indices = [[i - 1, j], [i, j - 1], [i, j + 1], [i + 1, j]]

    for indices in neighbour_indices:
        if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
            continue
        if original_grid.dynamic_recrystallization_number[indices[0]][indices[1]] == 1:
            recrystallized_grid[i][j] = 1
            # print(neighbours)
            return True
    return False

def cell_on_border(i, j):
    #m moore's neighbourhood
    neighbour_indices = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],
                         [i, j - 1],                     [i, j + 1],
                         [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]

    # neighbours = grid[i-1 : i+2, j-1 : j+2]
    
    # von-neuman neighbourhood
    #neighbour_indices = [[i - 1, j], [i, j - 1], [i, j + 1], [i + 1, j]]

    # if dynamic_recrystallization_number[i][j] == 1:
    #     border_grid[i][j] = 0
    for indices in neighbour_indices:
        if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
            continue
        # print('Cell Value, Neighbour Value : ', grid[i][j], grid[indices[0]][indices[1]])
        if original_grid.state[i][j] != original_grid.state[indices[0]][indices[1]]:
            # print(neighbours)
            border_grid[i][j] = 1
            return True
    return False

def main():
    global true_strain
    global true_stress
    global total_grains
    global original_grid
    global updation_grid

    original_grid.state = monte_carlo_initialization(grain_initialization())
    # original_grid.state = vornoi_initialization(original_grid.state)
    print(original_grid.state)
    # original_grid.dislocation_densities = initialize_dislocation_density()
    original_grid.orientation = initialize_orientation()
    
    k = 0
    N = states
    xDrx = []
    # setup the plot
    fig, ax = plt.subplots(1,1, figsize=(10,10))
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

    cb = plt.colorbar(update_grid_state(plt, fig, ax, cmap, bounds), spacing='proportional',ticks=bounds)
    cb.set_label('Custom cbar')
    

    update_grid_state(plt, fig, ax, cmap, bounds)
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

    while (k < iterations):
        prev_grid = copy.deepcopy(original_grid)
        counter_propagate1 = 0
        counter_propagate2 = 0
        i = int(random.randrange(n))
        j = int(random.randrange(n))
        # cell_on_border(i, j)
        # near_recrystallized_cell(i, j)
        if original_grid.dislocation_densities[i][j] >= P_cr:
            # if nucleatinon_step:
            #     if cell_on_border(i, j):
            #         ret = update_cell_state(i, j)
            #         update_indices += [i, j]
            #         nucleation_happened = True
            # else:
            #     if near_recrystallized_cell(i, j) and (original_grid.dynamic_recrystallization_number[i][j] == 0):
            #         ret = propagateGrainBoundary(i, j, k)
            #         propagate_indices += [i, j]
            #     else:
            #         updation_grid.dislocation_densities[i][j] = dislocation_energy(original_grid.dislocation_densities[i][j], original_grid.orientation[i][j])
            if original_grid.dynamic_recrystallization_number[i][j] == 0:
                if near_recrystallized_cell(i, j):
                    ret = propagateGrainBoundary(i, j, k)
                elif cell_on_border(i, j):
                    ret = update_cell_state(i, j)
                else:
                    updation_grid.dislocation_densities[i][j] = dislocation_energy(original_grid.dislocation_densities[i][j], original_grid.orientation[i][j])
            
                if ret != 1:
                    updation_grid.dislocation_densities[i][j] = dislocation_energy(original_grid.dislocation_densities[i][j], original_grid.orientation[i][j])
            else:
                updation_grid.dislocation_densities[i][j] = dislocation_energy(original_grid.dislocation_densities[i][j], original_grid.orientation[i][j])
        else:
            updation_grid.dislocation_densities[i][j] = dislocation_energy(original_grid.dislocation_densities[i][j], original_grid.orientation[i][j])

        for y in range(n):
            for x in range(n):
                if (y == i) and (j == x):
                    continue
                updation_grid.dislocation_densities[y][x] = dislocation_energy(original_grid.dislocation_densities[y][x], original_grid.orientation[y][x])

        cell_strain[i][j] = strain(original_grid.state[i][j])
        if nucleation_happened:
            nucleatinon_step = False
        original_grid = copy.deepcopy(updation_grid)
        # print(original_grid.state)
        # print("borer elements: ", counter_propagate1)
        # print("Not on border(after critical): ", counter_propagate2)

        unq, cnts = np.unique(original_grid.dynamic_recrystallization_number, return_counts=True)
        unq_array = dict(zip(unq, cnts))
        # print("updated cells: ", unq_array)
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
        N_c += (N_c - N_unchanged_grains)
        
        total_grains = N_c

        # print("Total Grains: ", total_grains)
        
        xDrx = xDrx + [float(N_r) / N_c]
        # print(xDrx)
        
        net_strain = np.mean(cell_strain)
        net_stress = simulated_stress()
        # print("Net strain value:", net_strain)
        # print("Net stress value:", net_stress)
        # true_strain = true_strain + [net_strain]
        true_strain += [strains[k]]
        true_stress += [net_stress]
        # print("True Stress: ", true_stress)
        print(k)
        # print(grid)
        # print(dislocation_densities)
        # pw.plot(true_strain, true_stress, pen='r')
        k += 1
        # print(grid)
        # BLOCK 1 START
        if k % 100 == 0:
            update_grid_state(plt, fig, ax, cmap, bounds)
            # ax[1].plot(true_strain, true_stress, 'b-')
            # ax[2].plot(true_strain, xDrx, 'r-')
            plt.pause(0.01)
        # try:
        #     if keyboard.is_pressed('q'):
        #         update_grid_state(plt, fig, ax, cmap, bounds)
        #         # ax[1].plot(true_strain, true_stress, 'b-')
        #         # ax[2].plot(true_strain, xDrx, 'r-')
        #         plt.pause(0.01)
        #     else:
        #         pass
        # except:
        #     pass
        # key = cv2.waitKey(5)
        
        # if 't' == chr(key & 255):
        #     print("Key Pressed")
        #     cv2.waitKey(0)  
        #     update_grid_state(plt, fig, ax, cmap, bounds)
        #     # ax[1].plot(true_strain, true_stress, 'b-')
        #     # ax[2].plot(true_strain, xDrx, 'r-')
        #     plt.pause(0.01)
        
        # BLOCK 1 END
    # BLOCK 2 START
    # update_grid_state(plt, fig, ax, cmap, bounds)
    # ax[1].plot(true_strain, true_stress, 'b-')
    # ax[2].plot(true_strain, xDrx, 'r-')
    #BLOCK 2 END
    plt.show()

if __name__ == "__main__":
    main()
