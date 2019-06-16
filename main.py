from helper import *

current_strain = 0
true_strain = []
true_stress = []
global new_grains
new_grains = 0

global total_grains
total_grains = n * n

global state_to_orientation
state_to_orientation = []

strains = interpolate_array(strains)
iterations = len(strains)
global k
k = 0


def main():
    original_grid = grid_class()
    updation_grid = copy.deepcopy(original_grid)

    

if __name__ == "__main__":
    main()