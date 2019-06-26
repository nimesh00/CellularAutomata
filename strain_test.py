import numpy as np


states = 50
n = 100
iterations = 30
strain_rate = 0.01
P_prev = 6.022 * (10 ** 14)
P_max = 8.214 * (10 ** 14)
P_cr = 8.214 * (10 ** 14)
critical_orientation = 15
k1 = 7.78 * 10 ** 8
k2 = 27.09
cd = 2 * 10 ** -6
# Mo = 1 * (10 ** 5)
Mo = 1.4 * (10 ** -6)
Tm = 1453  # degree kelvin
mu_o = 7.89 * 10 ** 4
Current_Temp = 1313
b = 2.49 * 10 ** -10
alpha = 0.5
orientation = 90
current_strain = 0.0172

def mu(temp):
    mu = mu_o * (1 - 0.64 * (temp - 27) / (Tm + 273))
    return mu

def mobility(orientation):
    return (Mo * (1 - np.exp( -5 * ((float(orientation) / critical_orientation) ** 4))))

def strain(orientation):
    # print('Cd: ', cd)
    # print('k1: ', k1)
    # print('k2: ', k2)
    # print('Mo: ', Mo)
    # print('orientation: ', orientation)
    # print('critical orientation: ', critical_orientation)
    delta_t = (cd * ((float(k2) / k1) ** 2)) / ((0.5 * mu(Current_Temp) * (b ** 2)) * mobility(orientation))
    # delta_t = (cd / ((0.5 * mu(Current_Temp) * (b ** 2)) * mobility(orientation)) * ((float(k2) / k1) ** 2))
    global current_strain
    delta_strain = strain_rate * delta_t
    current_strain +=  delta_strain
    print("delta strain: ", delta_strain)
    return current_strain, delta_strain

_, delta_strn = strain(orientation)
strain, _ = strain(orientation)

print(strain)

delta_t = (cd * ((float(k2) / k1) ** 2)) / ((1.28 * 10 ** -15) * mobility(orientation))
a = (0.01 * delta_t)
print(delta_strn)

delta_p = (k1 * np.sqrt(P_prev) - k2 * P_prev) * (delta_strn)
P_new = P_prev + delta_p
stress = alpha * mu(Current_Temp) * b * np.sqrt(P_new)
print(stress)