from globalVariables import *

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

    nucleation_probability = original_grid.dislocation_densities[i][j] * ((float(k2) / k1) ** 2)
    random_number = np.random.random()

    if random_number <= nucleation_probability:
        updation_grid.orientation[0][0][i][j] =                original_grid.orientation[0][0][favoured_indices[0]][favoured_indices[1]]
        updation_grid.orientation[0][1][i][j] =                original_grid.orientation[0][1][favoured_indices[0]][favoured_indices[1]]
        updation_grid.orientation[0][2][i][j] =                original_grid.orientation[0][2][favoured_indices[0]][favoured_indices[1]]
        updation_grid.state[i][j] =                            original_grid.state[favoured_indices[0]][favoured_indices[1]]
        updation_grid.dislocation_densities[i][j] =            P_RX_initial
        updation_grid.dynamic_recrystallization_number[i][j] = 1
        return 1
    return 0

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