import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt

#to develop the initial grain structure

iterate = 50
states = 10
n = 20
matx, maty = n, n
arr = np.zeros((matx,maty))
nx = []
ny = []
ns = []
for i in range(states):
    nx.append(random.randrange(matx))
    ny.append(random.randrange(maty))
    ns.append(int(random.randrange(states))+1)
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



mat = arr
#print(mat)

#to run the Monte Carlo simulation on the array obtained 
states_list = []

for i in range(states):
    states_list.append(i+1)

def most_common(lst):
    return max(set(lst), key=lst.count)

for p in range(iterate):
    for i in range(matx):
        for j in range(maty):
            not_similar_i = 0
            not_similar_f = 0
            switch = 0
            original = mat[i, j]
            final_decision = 0
            neighbourhood = []   #gives the position of the element in the list of states
            neighbour_indices = [[i - 1, j - 1], [i - 1, j], [i - 1, j + 1],
                                 [i, j - 1],                     [i, j + 1],
                                 [i + 1, j - 1], [i + 1, j], [i + 1, j + 1]]

            for indices in neighbour_indices:
                if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
                        continue
                neighbourhood.append(mat[indices[0],indices[1]])
                if mat[indices[0],indices[1]] == mat[i, j]:
                    continue
                else: 
                    not_similar_i += 1   #this would give me the total number of different cells in neighbour initially (or initial energy Ei)

                switch = most_common(neighbourhood)     #gives the element I would switch and calculate the energy with 
                
                   
            mat[i, j] = switch
            
            for indices in neighbour_indices:
                if (indices[0] < 0) or (indices[1] < 0) or (indices[0] > n - 1) or (indices[1] > n - 1):
                        continue
                neighbourhood.append(mat[indices[0],indices[1]])
                if mat[indices[0],indices[1]] == mat[i, j]:
                    continue
                else: 
                    not_similar_f += 1   #this would give me the total number of different cells in neighbour finally (or final energy Ef)


            #calculating whether switch must be accepted or not.

            delta_e = not_similar_f - not_similar_i

            if delta_e <= 0:
                final_decision = mat[i, j]
            else:
                prob = abs(delta_e)/8
                rand = random.random()
                if prob >= rand:
                    final_decision = mat[i,j]
                else:
                    final_decision = original

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
                    
                    























