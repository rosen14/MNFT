# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 16:37:40 2022

@author: Rosendo
"""

import numpy as np
import pandas as pd

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


def numerical_transitory_solution(n, phi_init, theta, t_max, deltaT, second_order_dirichlet = False):
    #Theta 0: Forwawrd Euler | 1/2: Crank Nicolson | 1: Backward Euler 
    
    
    h = L/n     # Tamaño de discretizacion (uniforme)

    x_caras = np.linspace(0, L, n + 1)  # Coordenada de caras
    
    # Coordenada de centroide de las celdas
    x_centroides = np.array([0.5*(x_caras[i+1] - x_caras[i]) + x_caras[i] \
                  for i in range(len(x_caras) - 1)])
    
    df_transitory = pd.DataFrame({'centroides': x_centroides,
                                  '0': phi_init})
    
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
        
    M_matrix = h/deltaT*np.identity(n) + K_matrix*theta
    Minv_matrix = np.linalg.inv(M_matrix)
    
    phi_k = phi_init
    t = 0
    while t < t_max:
        
        R_matrix = b_matrix + phi_k/deltaT*h - (1 - theta)*np.dot(K_matrix, phi_k)
    
        phi_k1 = np.dot(Minv_matrix, R_matrix)
        
        t = t + deltaT
        phi_k = phi_k1
        
        df_transitory[str(round(t,4))] = phi_k1
    
    return df_transitory


# Solucion de referencia ####
theta = 0.5
deltaT = 0.000001
t_max = 0.1

ref_solution = numerical_transitory_solution(n, phi_init, theta, t_max, deltaT, second_order_dirichlet = True)
###############################

dTBase = 0.001
discr_temp = [dTBase/2**(i-1) for i in range(1, 6)] 


error_df = pd.DataFrame()
for deltaT in discr_temp:
    phi_numerical_forward = numerical_transitory_solution(n, phi_init, 0, t_max, deltaT, second_order_dirichlet = True)
    phi_numerical_crank = numerical_transitory_solution(n, phi_init, 0.5, t_max, deltaT, second_order_dirichlet = True)
    phi_numerical_backward = numerical_transitory_solution(n, phi_init, 1, t_max, deltaT, second_order_dirichlet = True)
    
    norma1_forward = 1/n * np.sum(np.abs(ref_solution[str(t_max)] - phi_numerical_forward[str(t_max)]))
    norma1_crank = 1/n * np.sum(np.abs(ref_solution[str(t_max)] - phi_numerical_crank[str(t_max)]))
    norma1_backward = 1/n * np.sum(np.abs(ref_solution[str(t_max)] - phi_numerical_backward[str(t_max)]))
    
    error_df = error_df.append({
                    'discretizacion temporal': deltaT,
                    'norma 1 Forward': norma1_forward,
                    'norma 1 Crank': norma1_crank,
                    'norma 1 Backward': norma1_backward
                    },
                    ignore_index=True)
    
error_df['tasa error n1 Forward'] = pd.Series(np.log2(error_df['norma 1 Forward'][:-1].values/ \
                                                      error_df['norma 1 Forward'][1:].values))
error_df['tasa error n1 Crank'] = pd.Series(np.log2(error_df['norma 1 Crank'][:-1].values/ \
                                                      error_df['norma 1 Crank'][1:].values))
error_df['tasa error n1 Backward'] = pd.Series(np.log2(error_df['norma 1 Backward'][:-1].values/ \
                                                      error_df['norma 1 Backward'][1:].values))