from globalVariables import *

def cell_on_border(i, j):
    neighbour_indices = [                [i - 1, j],               
                         [i, j - 1],                     [i, j + 1],
                                         [i + 1, j],               ]

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
    neighbour_indices = [                [i - 1, j],               
                         [i, j - 1],                     [i, j + 1],
                                         [i + 1, j],               ]

    for indices in neighbour_indices:
        if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
            continue
        if original_grid.dynamic_recrystallization_number[indices[0]][indices[1]] == 1:
            recrystallized_grid[i][j] = 1
            return True
    return False

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