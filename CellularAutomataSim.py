'''
    COMMENT BLOCK 1 and UNCOMMENT BLOCK 2 : Updations first and then end plot. End results only
    UNCOMMENT BLOCK 1 and COMMENT BLOCK 2: Live updations after each iterations. Slower(can be optimized)
    Enter the correct values of the constants
'''
'''
Moore's neighbourhood used
Dislocation Density approach instead of Grain Boundary Velocity approach for determining nucleation

'''

import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import copy
import math

states = 10
n = 200
iterations = 30
strain_rate = 0.1
P_initial = 6.022 * (10 ** 14)
P_max = 8.214 * (10 ** 14)
P_cr = 8.214 * (10 ** 14)
critical_orientation = 7
k1 = 7.78 * 10 ** 8
k2 = 27.09
cd = 2
Mo = 1 * (10 ** 5)
Tm = 1453  # degree kelvin
mu_o = 7.89 * 10 ** 4
Current_Temp = 1313
b = 2.49 * 10 ** -10
alpha = 0.5

current_strain = 0
true_strain = []
true_stress = []
global new_grains
new_grains = 0

global total_grains
total_grains = n * n


def interpolate_array(array):
    #doubles the length
    current_length = len(array)
    print("original array: ", array)
    for i in range(current_length - 1):
        print(i, " elements added")
        print("New length: ", len(array))
        array.insert(2 * i + 1, ((array[2 * i] + array[2 * i + 1]) / 2))
        print("new array: ", array)
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
iterations = len(strains)
global k
k = 0

grid = np.random.randint(2, states, size = (n,n))
dislocation_energies = [[P_initial for l in range(n)] for m in range(n)]
dynamic_recrystallization_number = np.zeros((n, n))
distance_n = np.zeros((n, n))
cell_strain = np.zeros((n, n))


def mu(temp):
    mu = mu_o * (1 - 0.64 * (temp - 27) / (Tm + 273))
    return mu

def simulated_stress():
    # average_dislocation_energy = np.mean(dislocation_energies)
    global total_grains
    average_dislocation_energy = (np.mean(dislocation_energies) * n * n + P_initial * new_grains) / total_grains
    simulated_stress = alpha * mu(Current_Temp) * b * np.sqrt(average_dislocation_energy)
    return simulated_stress


def strain(orientation):
    # print('Cd: ', cd)
    # print('k1: ', k1)
    # print('k2: ', k2)
    # print('Mo: ', Mo)
    # print('orientation: ', orientation)
    # print('critical orientation: ', critical_orientation)
    delta_t = (cd * ((float(k2) / k1) ** 2)) / (0.5 * mu(Tm) * (b ** 2)) * (Mo * (1 - np.exp( -5 * ((float(orientation) / critical_orientation) ** 4))))
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

def update_state(P_prev, orientation_current):
    P_new = dislocation_energy(P_prev, orientation_current)
    new_state = orientation_current
    # print("Energy Difference", P_new - P_cr)
    global P_max
    P_max = np.max(dislocation_energies)
    # print("maximum DisDen: ", P_max)
    choice = 0
    if P_new < P_cr :
        new_state = orientation_current
    elif P_new >= P_cr:
        # print("P_new, P_max: ", P_new, P_max)
        nucleation_probability = P_new / P_max
        random_number = np.random.randint(0, 100)
        random_number = float(random_number) / 100
        # print("Random Number, Nucleation Prob", random_number, nucleation_probability)
        if random_number < nucleation_probability:
            choice = 1
            new_state = np.random.randint(2, states)

    return choice, new_state
    

def propagate_grain_boundary(y, x):
    neighbour_indices = [[y - 1, x - 1], [y - 1, x], [y - 1, x + 1],
                        [y , x - 1], [y , x], [y , x + 1],
                        [y + 1, x - 1], [y + 1, x], [y + 1, x + 1]]
    favoured_indices = neighbour_indices[0]
    min_delta_p = dislocation_energies[y][x] - dislocation_energies[favoured_indices[0]][favoured_indices[1]]
    for indices in neighbour_indices:
        if dynamic_recrystallization_number[indices[0]][indices[1]] == 1:
            delta_p = dislocation_energies[y][x] - dislocation_energies[indices[0]][indices[1]]
            if delta_p < min_delta_p:
                min_delta_p = delta_p
                favoured_indices = indices

    P_max = np.max(dislocation_energies)

    P_curr = dislocation_energies[y][x]

    nucleation_probability = P_curr / P_max

    random_number = np.random.randint(0, 99)
    random_number = float(random_number) / 100

    if random_number < nucleation_probability:
        grid[y][x] = grid[favoured_indices[0]][favoured_indices[1]]
        dislocation_energies[y][x] = 0
    else:
        dislocation_energies[y][x] = dislocation_energy(dislocation_energies[y][x], grid[y][x])


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
            tag[i*n + j] = grid[i][j]
    
    scat = ax.scatter(x,y,c=tag,cmap=cmap,     norm=norm)

    return scat

def check_mat_for_zero(matrix):
    height, width = matrix.shape

    for i in range(height):
        for j in range(width):
            if matrix[i][j] == 0:
                return True
    return False

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

def main():
    global true_strain
    global true_stress
    global grid
    global total_grains

    grid = vornoi_initialization(grid)

    
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
    
    while (k < iterations):
        prev_grid = copy.copy(grid)
        # number_border_elements = 0
        for i in range(n):
            for j in range(n):

                if dynamic_recrystallization_number[i][j] == 1:
                    continue

                # Moore's neighbourhood
                neighbour_indices = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],
                                     [i , j - 1], [i , j], [i , j + 1],
                                     [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]
                # print(neighbour_indices)
                border = False
                near_nucleus = False
                for indices in neighbour_indices:
                    if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
                        continue
                    # print('Cell Value, Neighbour Value : ', grid[i][j], grid[indices[0]][indices[1]])
                    if grid[i][j] != grid[indices[0]][indices[1]]:
                        border = True
                    if dynamic_recrystallization_number[indices[0]][indices[1]] == 1:
                        near_nucleus = True
                # if border is False:
                #     number_border_elements += 1
                #     print(number_border_elements)
                if border is True:
                    if near_nucleus is True:
                        propagate_grain_boundary(i, j)
                    else:
                        ret, grid[i][j] = update_state(dislocation_energies[i][j], grid[i][j])
                        if ret == 1:
                            dislocation_energies[i][j] = 0
                            dynamic_recrystallization_number[i][j] = 1
                            distance_n[i][j] = 0
                        else:
                            dislocation_energies[i][j] = dislocation_energy(dislocation_energies[i][j], grid[i][j])
                else:
                    dislocation_energies[i][j] = dislocation_energy(dislocation_energies[i][j], grid[i][j])
                cell_strain[i][j] = strain(grid[i][j])
        
        unq, cnts = np.unique(dislocation_energies, return_counts=True)
        unq_array = dict(zip(unq, cnts))

        print("Dslctn dnsts: ", unq_array)

        difference_matrix = grid - prev_grid
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

        print("Total Grains: ", total_grains)
        
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
        print("Iteration",k, "completed")
        # print(grid)
        # print(dislocation_energies)
        # pw.plot(true_strain, true_stress, pen='r')
        k += 1
        # print(grid)
        # BLOCK 1 START
        
        # update_grid_state(plt, fig, ax[0], cmap, bounds)
        # ax[1].plot(true_strain, true_stress, 'b-')
        # ax[2].plot(true_strain, xDrx, 'r-')
        # plt.pause(0.01)
        
        # BLOCK 1 END
    # BLOCK 2 START
    update_grid_state(plt, fig, ax[0], cmap, bounds)
    ax[1].plot(true_strain, true_stress, 'b-')
    ax[2].plot(true_strain, xDrx, 'r-')
    #BLOCK 2 END
    plt.show()

if __name__ == "__main__":
    main()
