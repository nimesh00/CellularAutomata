from globalVariables import *
from helperFunctions import *
from workingFunctions import *
from initializationFunctions import *
from constants import *

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

    while (k < iterations):
        if np.sum(original_grid.dynamic_recrystallization_number) == n * n:
            break
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
                            misorientations += [mis_orientation(original_grid.orientation[0][0][i][j], original_grid.orientation[0][1][i][j], original_grid.orientation[0][2][i][j], original_grid.orientation[0][0][indices[0]][indices[1]], original_grid.orientation[0][1][indices[0]][indices[1]], original_grid.orientation[0][2][indices[0]][indices[1]])]
                            delta_dislocations += [abs(original_grid.dislocation_densities[i][j] - original_grid.dislocation_densities[indices[0]][indices[1]])]

        for i in range(len(misorientations)):
            velocities += [grain_velocity(delta_dislocations[i], misorientations[i])]
            
        _, delta_strain = strain(np.mean(velocities))


        #your normal code continues from here
        for i in range(n):
            for j in range(n):
                cell_on_border(i, j)
                near_recrystallized_cell(i, j)
                if original_grid.dislocation_densities[i][j] >= P_cr:
                    if nucleatinon_step:
                        if cell_on_border(i, j):
                            ret = update_recrystallized_cell_state(i, j)                                                                
                            update_indices += [i, j]
                            nucleation_happened = True
                    else:
                        if near_recrystallized_cell(i, j) and (original_grid.dynamic_recrystallization_number[i][j] == 0):
                            ret = propagateGrainBoundary(i, j, k)
                            propagate_indices += [i, j]
                        else:
                            updation_grid.dislocation_densities[i][j] = update_dislocation_density(updation_grid.dislocation_densities[i][j], delta_strain)
                    if ret != 1:
                        updation_grid.dislocation_densities[i][j] = update_dislocation_density(updation_grid.dislocation_densities[i][j], delta_strain)
                else:
                    updation_grid.dislocation_densities[i][j] = update_dislocation_density(updation_grid.dislocation_densities[i][j], delta_strain)
                
                # if original_grid.dislocation_densities[i][j] >= P_cr:
                #     if original_grid.dynamic_recrystallization_number[i][j] == 0:
                #         if near_recrystallized_cell(i, j):
                #             ret = propagateGrainBoundary(i, j, k)
                #         elif cell_on_border(i, j):
                #             ret = update_recrystallized_cell_state(i, j)
                #         else:
                #             updation_grid.dislocation_densities[i][j] = update_dislocation_density(updation_grid.dislocation_densities[i][j], delta_strain)
                #         if ret != 1:
                #             updation_grid.dislocation_densities[i][j] = update_dislocation_density(updation_grid.dislocation_densities[i][j], delta_strain)
                #     else:
                #         updation_grid.dislocation_densities[i][j] = update_dislocation_density(updation_grid.dislocation_densities[i][j], delta_strain)
                # else:
                #     updation_grid.dislocation_densities[i][j] = update_dislocation_density(updation_grid.dislocation_densities[i][j], delta_strain)
                #cell_strain[i][j], _ = strain(original_grid.state[i][j])
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