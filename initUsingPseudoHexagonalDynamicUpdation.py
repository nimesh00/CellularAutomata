import copy
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

states = 10
n = 100

class grid_class:
    def __init__(self):
        self.state = np.random.randint(2, states, size = (n,n))
        self.dynamic_recrystallization_number = np.zeros((n, n))

original_grid = grid_class()
original_grid.state = np.zeros((n, n))
updation_grid = copy.deepcopy(original_grid)

def pseudoHexagonalNeighbourInit():
    no_of_nucleus = np.random.randint(states, 2 * states)
    for k in range(no_of_nucleus):
        Y = np.random.randint(0, n)
        X = np.random.randint(0, n)
        S = np.random.randint(2, states)
        updation_grid.state[Y][X] = S
        updation_grid.dynamic_recrystallization_number[Y][X] = 1

    original_grid = copy.deepcopy(updation_grid)

    flag = 0
    m = 0
    while np.sum(original_grid.dynamic_recrystallization_number) != (n * n):
        for i in range(n):
            for j in range(n):
                near_recrystallized_cell = False
                if flag == 0:
                    neighbour_indices = [[i - 1, j - 1], [i - 1, j], 
                                        [i, j - 1],                      [i, j + 1],
                                                        [i + 1, j], [i + 1, j + 1]]
                    flag = 1
                elif flag == 1:
                    neighbour_indices = [                [i - 1, j], [i - 1, j + 1],
                                        [i, j - 1],                      [i, j + 1],
                                        [i + 1, j - 1], [i + 1, j],               ]
                    flag = 0
                
                for indices in neighbour_indices:
                    if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
                        continue
                    if original_grid.dynamic_recrystallization_number[indices[0]][indices[1]] == 1:
                        near_recrystallized_cell = True
                        updation_grid.state[i][j] = original_grid.state[indices[0]][indices[1]]
                        updation_grid.dynamic_recrystallization_number[i][j] = 1
                        break
        original_grid = copy.deepcopy(updation_grid)
        print("Iteration", m)
        m += 1
    return original_grid.state


def main():

    no_of_nucleus = np.random.randint(states, 2 * states)
    for k in range(no_of_nucleus):
        Y = np.random.randint(0, n)
        X = np.random.randint(0, n)
        S = np.random.randint(2, states)
        updation_grid.state[Y][X] = S
        updation_grid.dynamic_recrystallization_number[Y][X] = 1

    original_grid = copy.deepcopy(updation_grid)

    flag = 0
    m = 0
    mat = np.zeros((n, n))
    N = states
    fig, ax = plt.subplots(1,1, figsize=(10,10))

    # define the colormap
    cmap = plt.cm.jet
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

    # define the bins and normalize
    bounds = np.linspace(0,N,N+1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    x = np.zeros((n*n))
    y = np.zeros((n*n))
    tag = np.zeros(n*n)
    for i in range(n):
        for j in range(n):
            x[i*n + j] = i
            y[i*n + j] = j
            tag[i*n + j] = mat[i][j]

    # make the scatter
    scat = ax.scatter(x, y, marker='s', c=tag,cmap=cmap,     norm=norm)
    cb = plt.colorbar(scat, spacing='proportional',ticks=bounds)
    cb.set_label('Custom cbar')

    while np.sum(original_grid.dynamic_recrystallization_number) != (n * n):
        for i in range(n):
            for j in range(n):
                if original_grid.dynamic_recrystallization_number[i][j] == 1:
                    continue
                near_recrystallized_cell = False
                if flag == 0:
                    neighbour_indices = [[i - 1, j - 1], [i - 1, j], 
                                        [i, j - 1],                      [i, j + 1],
                                                        [i + 1, j], [i + 1, j + 1]]
                    # flag = 1
                elif flag == 1:
                    neighbour_indices = [                [i - 1, j], [i - 1, j + 1],
                                        [i, j - 1],                      [i, j + 1],
                                        [i + 1, j - 1], [i + 1, j],               ]
                    # flag = 0
                
                for indices in neighbour_indices:
                    if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
                        continue
                    if original_grid.dynamic_recrystallization_number[indices[0]][indices[1]] == 1:
                        near_recrystallized_cell = True
                        updation_grid.state[i][j] = original_grid.state[indices[0]][indices[1]]
                        updation_grid.dynamic_recrystallization_number[i][j] = 1
                        break
        flag = not flag
        original_grid = copy.deepcopy(updation_grid)
        mat = original_grid.state
        print("Iteration", m)
        m += 1
            # define the data
        x = np.zeros((n*n))
        y = np.zeros((n*n))
        tag = np.zeros(n*n)
        for i in range(n):
            for j in range(n):
                x[i*n + j] = i
                y[i*n + j] = j
                tag[i*n + j] = mat[i][j]

        # make the scatter
        scat = ax.scatter(x, y, marker='s', c=tag,cmap=cmap,     norm=norm)
        plt.pause(0.01)

    
    





if __name__ == "__main__":
    main()