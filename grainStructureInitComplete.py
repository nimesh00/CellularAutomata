import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt

#to develop the initial grain structure

iterations = 100000
states = 50
n = 100

mat = np.zeros((n,n))
for m in range(states):
    i = int(random.randrange(n))
    j = int(random.randrange(n))
    mat[i][j] = 1 + random.randrange(states)

#print(mat)

for m in range(iterations):
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


#photo of the thing
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
scat = ax.scatter(x,y, marker='s', c=tag,cmap=cmap,     norm=norm)
# create the colorbar
cb = plt.colorbar(scat, spacing='proportional',ticks=bounds)
cb.set_label('Custom cbar')
plt.pause(2) 


#to run the Monte Carlo simulation on the array obtained 
states_list = []

for i in range(states):
    states_list.append(i+1)

'''def most_common(lst):
    return max(set(lst), key=lst.count)'''

def most_common(List): 
    counter = 0
    num = List[0] 
      
    for i in List: 
        curr_frequency = List.count(i) 
        if(curr_frequency> counter): 
            counter = curr_frequency 
            num = i 
  
    return num

for m in range(iterations):
    print(m)
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



#photo of the thing
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
scat = ax.scatter(x,y, marker='s', c=tag,cmap=cmap,     norm=norm)
# create the colorbar
cb = plt.colorbar(scat, spacing='proportional',ticks=bounds)
cb.set_label('Custom cbar')
plt.show() 

            
