import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import copy
import math


for i in range(100):
    ph = np.random.randint(360, size=(2, 1)) 
    th = np.random.randint(180, size=(2, 1))
    omega = np.random.randint(360, size=(2, 1))

    R_1 = [[math.cos(omega[0][0]) * math.cos(ph[0][0]) - math.sin(omega[0][0]) * math.sin(ph[0][0]) * math.cos(th[0][0]),    math.cos(omega[0][0]) * math.sin(ph[0][0]) + math.sin(omega[0][0]) * math.cos(ph[0][0]) * math.cos(th[0][0]),   math.sin(omega[0][0]) * math.sin(th[0][0])],
           [-math.sin(omega[0][0]) * math.cos(ph[0][0]) - math.cos(omega[0][0]) * math.sin(ph[0][0]) * math.cos(th[0][0]),  -math.sin(ph[0][0]) * math.sin(omega[0][0]) + math.cos(ph[0][0]) * math.cos(omega[0][0]) * math.cos(th[0][0]),   math.cos(omega[0][0]) * math.sin(th[0][0])],
           [math.sin(ph[0][0]) * math.sin(th[0][0]),                                                                             math.sin(th[0][0]) * math.cos(ph[0][0]),                                                                                                math.cos(th[0][0])]]

    R_2 = [[math.cos(omega[1][0]) * math.cos(ph[1][0]) - math.sin(omega[1][0]) * math.sin(ph[1][0]) * math.cos(th[1][0]),   math.cos(omega[1][0]) * math.sin(ph[1][0]) +  math.sin(omega[1][0]) * math.cos(ph[1][0]) * math.cos(th[1][0]),  math.sin(omega[1][0]) * math.sin(th[1][0])],
           [-math.sin(omega[1][0]) * math.cos(ph[1][0]) - math.cos(omega[1][0]) * math.sin(ph[1][0]) * math.cos(th[1][0]),  -math.sin(ph[1][0]) * math.sin(omega[1][0]) + math.cos(ph[1][0]) * math.cos(omega[1][0]) * math.cos(th[1][0]),  math.cos(omega[1][0]) * math.sin(th[1][0])],
           [math.sin(ph[1][0]) * math.sin(th[1][0]),                                                                          math.sin(th[1][0]) * math.cos(ph[1][0]),                                                                                                math.cos(th[1])]]

    #print(R_2) 

    trans_R_2 = [[R_2[j][i] for j in range(len(R_2))] for i in range(len(R_2[0]))]

    #print(trans_R_2) 

    mul = np.zeros([3, 3])

    for i in range(3):
        for j in range(3):
            mul[i][j] = R_1[i][j] * trans_R_2[i][j]

    #print(mul) 

    trace = 0

    for i in range(3):
        for j in range(3):
            if i == j:
                trace += mul[i][j]

    final = (trace - 1) / 2
    #print(final)

    inv = np.arccos(final) * 180 / math.pi
    print(inv)
    
