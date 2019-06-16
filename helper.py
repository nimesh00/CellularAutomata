import numpy as np

states = 20
n = 100
iterations = 30
strain_rate = 0.01
P_initial = 6.022 * (10 ** 14)
P_max = 8.214 * (10 ** 14)
P_cr = 7.214 * (10 ** 14)
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

class grid_class:
    def __init__(self):
        self.state = np.random.randint(2, states, size = (n,n))
        self.dislocation_densities = [[P_initial for l in range(n)] for m in range(n)]
        self.dynamic_recrystallization_number = np.zeros((n, n))
        self.orientation = np.random.randint(1, 180, size = (n,n))

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
        
def monte_carlo_initialization(mat):
    iterations = 100000
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

def cell_on_border(i, j):
    #m moore's neighbourhood
    neighbour_indices = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],
                        [i, j - 1],                      [i, j + 1],
                        [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]

    # neighbours = grid[i-1 : i+2, j-1 : j+2]
    
    # von-neuman neighbourhood
    # neighbour_indices = [[i - 1, j], [i, j - 1], [i, j + 1], [i + 1, j]]

    # if dynamic_recrystallization_number[i][j] == 1:
    #     border_grid[i][j] = 0
    for indices in neighbour_indices:
        if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
            continue
        # print('Cell Value, Neighbour Value : ', grid[i][j], grid[indices[0]][indices[1]])
        if grid[i][j] != grid[indices[0]][indices[1]]:
            # print(neighbours)
            border_grid[i][j] = 1
            return True
    return False
        
def near_recrystallized_cell(i, j):
    #m moore's neighbourhood
    neighbour_indices = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],
                        [i, j - 1],                      [i, j + 1],
                        [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]
    
    # neighbours = dynamic_recrystallization_number[i-1 : i+2, j-1 : j+2]
    
    # von-neuman neighbourhood
    # neighbour_indices = [[i - 1, j], [i, j - 1], [i, j + 1], [i + 1, j]]

    for indices in neighbour_indices:
        if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
            continue
        if dynamic_recrystallization_number[indices[0]][indices[1]] == 1:
            recrystallized_grid[i][j] = 1
            # print(neighbours)
            return True
    return False
        
def mu(temp):
    mu = mu_o * (1 - 0.64 * (temp - 27) / (Tm + 273))
    return mu

def simulated_stress():
    # average_dislocation_energy = np.mean(dislocation_densities)
    global total_grains
    average_dislocation_energy = np.mean(dislocation_densities)
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
