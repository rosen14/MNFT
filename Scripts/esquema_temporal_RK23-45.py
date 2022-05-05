# -*- coding: utf-8 -*-
"""
Created on Wed May  4 21:13:05 2022

@author: Rosendo
"""

import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp

# CONSTANTES ----------------------------------------------------------------#
L = 1        # Largo de barra [m]
S = 1        # Área transversal [m2]
Q = 1     # Fuente/(rho*Cp) [°C/s]
alpha = 1 # Difusividad [m2/s]  k/(rho*Cp)
n = 10       # 

# Condiciones de Borde de Dirichlet (sobre la temperatura) ------------------#
phi_0 = 0
phi_L = 0

# Condicion inicial ---------------------------------------------------------#
phi_init = np.zeros(n)


def get_matrix(n, second_order_dirichlet = False):

    h = L/n     # Tamaño de discretizacion (uniforme)
    
    b_matrix = np.zeros((n,))  #Matriz término fuente y C.B.
    K_matrix = np.zeros((n,n)) #Matriz de difusión
    
    
    # Por condición de borde
    if second_order_dirichlet:
        K_matrix[0, 0:2] = -alpha/h*np.array([-4, 4/3])
        K_matrix[n-1, n-2:n] = -alpha/h*np.array([4/3, -4])
        
        b_matrix[0] = Q*h + 8/3*alpha/h*phi_0
        b_matrix[n-1] = Q*h + 8/3*alpha/h*phi_L
    else:
        K_matrix[0, 0:2] = -alpha/h*np.array([-3, 1])
        K_matrix[n-1, n-2:n] = -alpha/h*np.array([1, -3])
        
        b_matrix[0] = Q*h + 2*alpha/h*phi_0
        b_matrix[n-1] = Q*h + 2*alpha/h*phi_L

    for i in range(1, n-1):
        K_matrix[i,i-1:i+2] = -alpha/h*np.array([1, -2, 1])
        b_matrix[i] = Q*h

    return K_matrix, b_matrix

K_matrix, b_matrix = get_matrix(n, second_order_dirichlet = True)

def f(t, y):
    return np.dot(-K_matrix, y) + b_matrix

t0 = 0
tf = 0.1 * 10
#Pareciera haber un bug de SciPy y tengo que multiplicar el tmax por 10
sol_K23 = solve_ivp(f, [t0, tf], phi_init, 'RK23')
sol_K45 = solve_ivp(f, [t0, tf], phi_init, 'RK45')

print('Nro. pasos temporales: ', sol_K23.t.size)

ref_solution = np.array([0.01626624, 0.04184973, 0.05958297, 0.07066941, 0.07598301,
       0.07598301, 0.07066941, 0.05958297, 0.04184973, 0.01626624])

df_comp = pd.DataFrame({'REF': ref_solution,
                   'K23 err': ref_solution - sol_K23.y[:, -1],
                   'K45 err': ref_solution - sol_K45.y[:, -1]})