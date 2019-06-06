import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv
import pandas as pd

states = 20
n = 920
iterations = 1000000

# mat = np.zeros((n,n))
# for i in range(n):
#     for j in range(n):
#         mat[i][j] = random.randint(1,states)
# mat

# for i in range(iterations):
#     x = random.randint(0,n) - 1
#     y = random.randint(0,n) - 1
#     arr = [0]*states
#     if(x>0 and y>0):
#         arr[int(mat[x-1][y-1]) - 1] += 1
#     if(x>0):
#         arr[int(mat[x-1][y]) - 1] += 1
#     if(x>0 and y<n-1):
#         arr[int(mat[x-1][y+1]) - 1] += 1
#     if(y>0):
#         arr[int(mat[x][y-1]) - 1] += 1
#     if(y<n-1):
#         arr[int(mat[x][y+1]) - 1] += 1
#     if(x<n-1 and y>0):
#         arr[int(mat[x+1][y-1]) - 1] += 1
#     if(x<n-1):
#         arr[int(mat[x+1][y]) - 1] += 1
#     if(x<n-1 and y<n-1):
#         arr[int(mat[x+1][y+1]) - 1] += 1
#     temp_max = -10000
#     index = 0
#     for j in range(states):
#         if(arr[j]>temp_max):
#             temp_max = arr[j]
#             index = j+1
#     mat[x][y] = index

matrix_file = open('monteCarlo_init_matrix.csv', 'r')
data = matrix_file.read()
print(matrix_file)
# for i in range(len(matrix_file)):
#         mat[i] = int(matrix_file[0])
# mat = pd.read_csv(r'monteCarlo_init_matrix.csv')        
        
print(mat)

N = states

# setup the plot
fig, ax = plt.subplots(1,1, figsize=(10,10))
# define the data
x = np.zeros((n*n))
y = np.zeros((n*n))
tag = np.zeros(n*n)
for i in range(n):
    for j in range(n):
        x[i*n + j] = i
        y[i*n + j] = j
        tag[i*n + j] = mat[i][j]


# define the colormap
cmap = plt.cm.jet
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

# define the bins and normalize
bounds = np.linspace(0,N,N+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

# make the scatter
scat = ax.scatter(x,y,c=tag,cmap=cmap,     norm=norm)
# create the colorbar
cb = plt.colorbar(scat, spacing='proportional',ticks=bounds)
cb.set_label('Custom cbar')
plt.show()