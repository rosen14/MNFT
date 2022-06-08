# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 15:56:50 2022

@author: Rosendo
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot
from plotly.tools import make_subplots

# CONSTANTES ----------------------------------------------------------------#
L = 1        # Largo de barra [m]
S = 1        # Área transversal [m2]
Q = 0    # Fuente/(rho*Cp) [°C/s]
alpha = 0 # Difusividad [m2/s]  k/(rho*Cp)
n = 20       # 
t_max = 0.5    # Tiempo máximo en [s]
V = 1
# Condiciones de Borde de Dirichlet (sobre la temperatura) ------------------#
phi_0 = 1
phi_L = 2

# Condicion inicial ---------------------------------------------------------#
phi_init = np.zeros(n)

#Theta 0: Forwawrd Euler | 1/2: Crank Nicolson | 1: Backward Euler 

def analytic_solution(x):
    c0 = phi_0
    c1 = phi_L/L + Q*L/(2*alpha) - c0/L
    return -Q*np.power(x,2)/(2*alpha) + c1*x + c0

def numerical_transitory_solution(n, phi_init, theta, t_max, deltaT, upwind = True, second_order_dirichlet = False):
    h = L/n     # Tamaño de discretizacion (uniforme)

    x_caras = np.linspace(0, L, n + 1)  # Coordenada de caras
    
    # Coordenada de centroide de las celdas
    x_centroides = np.array([0.5*(x_caras[i+1] - x_caras[i]) + x_caras[i] \
                  for i in range(len(x_caras) - 1)])
    
    df_transitory = pd.DataFrame({'centroides': x_centroides,
                                  '0': phi_init})
    
    b_matrix = np.zeros((n,))  #Matriz término fuente y C.B.
    K_matrix = np.zeros((n,n)) #Matriz de difusión
    
    ba_matrix = np.zeros((n,))
    A_matrix = np.zeros((n,n)) #Matriz de advección
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
        
    if upwind:
        A_matrix[0, 0] = V
        ba_matrix[0] = V*phi_0
        for i in range(1, n):
            A_matrix[i,i-1:i+1] = [-V, V]
    
    else: #central differences
        A_matrix[0, 0:2] = [V/2, V/2]
        A_matrix[n-1, n-2:n] = [-V/2, -V/2]
        ba_matrix[0] = V*phi_0
        ba_matrix[n-1] = -V*phi_L
        for i in range(1, n-1):
            A_matrix[i,i-1:i+2] = [-V/2, 0, V/2]
        
    K_matrix += A_matrix
    b_matrix += ba_matrix
    M_matrix = h/deltaT*np.identity(n) + K_matrix*theta
    Minv_matrix = np.linalg.inv(M_matrix)
    
    phi_k = phi_init
    t = 0
    precision = int(-np.log10(deltaT))
    while t < t_max:
        
        R_matrix = b_matrix + phi_k/deltaT*h - (1 - theta)*np.dot(K_matrix, phi_k)
    
        phi_k1 = np.dot(Minv_matrix, R_matrix)
        
        t = round(t + deltaT, precision)
        phi_k = phi_k1
        
        df_transitory[str(t)] = phi_k1
    
    return df_transitory

t_max = 0.5
df_transitory = numerical_transitory_solution(n, phi_init, 0.5, t_max, 0.01, upwind = False, second_order_dirichlet = False)
df_transitory[str(t_max)]