from globalVariables import *

def initialize_dislocation_density():
    possible_random_numbers = np.random.random((states+1))
    dislocation_densities = [[P_initial for i in range(n)] for j in range(n)]                                       #creates states-2 diff random values only?
    for i in range(n):
        for j in range(n):
            # dislocation_densities[i][j] = (2 * possible_random_numbers[original_grid.state[i][j]] + 6) * (10 ** 14)
            dislocation_densities[i][j] = (2 * np.random.random() + 5) * (10 ** 14)
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