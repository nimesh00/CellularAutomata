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
cell_strain = np.zeros((n, n))

def mu(temp):
    mu = mu_o * (1 - 0.64 * (temp - 27) / (Tm + 273))
    return mu

def simulated_stress():
    average_dislocation_energy = np.mean(dislocation_energies)
    simulated_stress = alpha * mu(Current_Temp) * b * np.sqrt(average_dislocation_energy)
    return simulated_stress


def strain(orientation):
    # print('Cd: ', cd)
    # print('k1: ', k1)
    # print('k2: ', k2)
    # print('Mo: ', Mo)
    # print('orientation: ', orientation)
    # print('critical orientation: ', critical_orientation)
    delta_t = (cd * ((float(k2) / k1) ** 2)) / (Mo * (1 - np.exp( -5 * ((float(orientation) / critical_orientation) ** 4))))
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
    if P_new < P_cr :
        new_state = orientation_current
    elif P_new >= P_cr:
        nucleation_probability = P_new / P_max
        # choice = np.random.choice(2, 1, p = [1 - nucleation_probability, nucleation_probability])
        P_nucleate = int(10.0 * nucleation_probability)
        # print('P_nucleate: ', P_nucleate)
        choice_array = [1 for i in range(P_nucleate)]
        # num_zeroes = 10 - len(choice_array)
        for j in range(10):
            if (j + 1 )> P_nucleate:
                choice_array += [0]
        # print(choice_array)
        # choice_array = [0, 1, 0, 0, 1, 1, 0, 1, 1, 1]
        random_index = np.random.randint(0, 10)
        # print(random_index)
        choice = choice_array[random_index]
        if choice == 1 :
            new_state = np.random.randint(2, states)

    return new_state
    

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


# def updating_plot(pw, x, y):
#     pw = pg.plot(x, y, pen = 'r')

def check_mat_for_zero(matrix):
    height, width = matrix.shape

    for i in range(height):
        for j in range(width):
            if matrix[i][j] == 0:
                return True
    return False

def monte_carlo_initialization(mat):
    iterations = 700000
    for i in range(iterations):
        x = random.randint(0,n) - 1
        y = random.randint(0,n) - 1
        arr = [0]*states
        if(x>0 and y>0):
            arr[int(mat[x-1][y-1]) - 1] += 1
        if(x>0):
            arr[int(mat[x-1][y]) - 1] += 1
        if(x>0 and y<n-1):
            arr[int(mat[x-1][y+1]) - 1] += 1
        if(y>0):
            arr[int(mat[x][y-1]) - 1] += 1
        if(y<n-1):
            arr[int(mat[x][y+1]) - 1] += 1
        if(x<n-1 and y>0):
            arr[int(mat[x+1][y-1]) - 1] += 1
        if(x<n-1):
            arr[int(mat[x+1][y]) - 1] += 1
        if(x<n-1 and y<n-1):
            arr[int(mat[x+1][y+1]) - 1] += 1
        temp_max = -10000
        index = 1
        # print(arr)
        for j in range(states):
            if(arr[j]>temp_max):
                temp_max = arr[j]
                index = j+1
        # print(index)
        mat[x][y] = index

    # print(check_mat_for_zero(mat))
    # print(mat)
        
    return mat

def vornoi_initialization(mat):
    import math
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
    # grid = monte_carlo_initialization(grid)
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
    
    while (k < iterations):
        prev_grid = copy.copy(grid)
        number_border_elements = 0
        for i in range(n):
            for j in range(n):
                neighbour_indices = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],
                                     [i , j - 1], [i , j], [i , j + 1],
                                     [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]
                # print(neighbour_indices)
                border = False
                for indices in neighbour_indices:
                    if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
                        continue
                    # print('Cell Value, Neighbour Value : ', grid[i][j], grid[indices[0]][indices[1]])
                    if grid[i][j] != grid[indices[0]][indices[1]]:
                        border = True
                # if border is False:
                #     number_border_elements += 1
                #     print(number_border_elements)
                if border is True:
                    grid[i][j] = update_state(dislocation_energies[i][j], grid[i][j])
                    dislocation_energies[i][j] = dislocation_energy(dislocation_energies[i][j], grid[i][j])
                    cell_strain[i][j] = strain(grid[i][j])
        
        difference_matrix = grid - prev_grid
        # N_unchanged_grains = difference_matrix.count(0)
        unique, counts = np.unique(difference_matrix, return_counts=True)
        unique_array = dict(zip(unique, counts))
        print(unique_array)
        N_unchanged_grains = unique_array[0]

        # total number of grains
        N_c =  n * n
        N_r = N_c - N_unchanged_grains
        
        xDrx = xDrx + [float(N_r) / N_c]
        print(xDrx)
        
        net_strain = np.mean(cell_strain)
        net_stress = simulated_stress()
        # print("Net strain value:", net_strain)
        # print("Net stress value:", net_stress)
        true_strain += [strains[k]]
        true_stress = true_stress + [net_stress]
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
