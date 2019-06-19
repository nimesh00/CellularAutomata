import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt

iterate = 100
states = 20
n = 50
N = states
matx, maty = n, n
arr = np.zeros((matx,maty))
nx = []
ny = []
ns = []
for i in range(states):
    nx.append(random.randrange(matx))
    ny.append(random.randrange(maty))
    ns.append(int(random.randrange(states+1)))
    arr[nx[i],ny[i]] = ns[i]

#print(arr)
#print(ns)        
for p in range(iterate):
    for i in range(matx):
        for j in range(maty):
            if p%2 == 0:
                neighbour_indices = [[i - 1, j - 1], [i - 1, j], 
                                    [i, j - 1],                      [i, j + 1],
                                                    [i + 1, j], [i + 1, j + 1]]

            else:
                neighbour_indices = [                [i - 1, j], [i - 1, j + 1],
                                    [i, j - 1],                      [i, j + 1],
                                    [i + 1, j - 1], [i + 1, j],               ]

            if arr[i,j] in ns:
                for indices in neighbour_indices:
                    if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
                        continue
                    if arr[indices[0]][indices[1]] in ns:
                        continue
                    else:
                        arr[indices[0]][indices[1]] = arr[i,j]


print(arr)

mat = arr

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

