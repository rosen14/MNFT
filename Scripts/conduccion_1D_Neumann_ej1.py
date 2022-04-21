# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 20:56:33 2022

@author: Rosendo
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot

# CONSTANTES ----------------------------------------------------------------#
L = 10        # Largo de barra [m]
S = 1        # Área transversal [m2]
Q = 1     # Fuente/(rho*Cp) [°C/s]
alpha = 1 # Difusividad [m2/s]  k/(rho*Cp)


# Condiciones de Borde de Dirichlet (sobre la temperatura) ------------------#
phi_0 = 2
grad_L = -2

def analytic_solution(x):
    c0 = phi_0
    c1 = grad_L+Q*L/alpha
    return -Q*np.power(x,2)/(2*alpha) + c1*x + c0


def numerical_solution(n):
    h = L/n     # Tamaño de discretizacion (uniforme)
    
    x_caras = np.linspace(0, L, n + 1)  # Coordenada de caras
    
    # Coordenada de centroide de las celdas
    x_centroides = np.array([0.5*(x_caras[i+1] - x_caras[i]) + x_caras[i] \
                  for i in range(len(x_caras) - 1)])
    
    b_matrix = np.zeros((n,))  #Matriz término fuente y C.B.
    K_matrix = np.zeros((n,n)) #Matriz de difusión
    
    # Por condición de borde
    K_matrix[0, 0:2] = -alpha/h*np.array([-3, 1])
    K_matrix[n-1, n-2:n] = -alpha/h*np.array([1, -1])
    
    b_matrix[0] = Q*h + 2*alpha/h*phi_0
    b_matrix[n-1] = Q*h + alpha*grad_L

    for i in range(1, n-1):
        K_matrix[i,i-1:i+2] = -alpha/h*np.array([1, -2, 1])
        b_matrix[i] = Q*h
        
    phi_numerical = np.linalg.solve(K_matrix, b_matrix)

    return x_centroides, phi_numerical

# PLOT ----------------------------------------------------------------------#
fig = go.Figure()

xc, phi_numeric = numerical_solution(n = 10)
x_domain = xc
phi_analytic = analytic_solution(np.concatenate(([0], xc, [L])))

# Add traces
fig.add_trace(
    go.Scatter(
        x = x_domain,
        y = phi_numeric,
        mode='lines+markers', name = 'Numerical',
        marker=dict(
        size=8))
    )
fig.add_trace(
    go.Scatter(
        x = np.concatenate(([0], xc, [L])),
        y = phi_analytic,
        mode='lines+markers', name = 'Analytic'))

fig.update_layout(title='Difusión estacionario 1D - Problema condición Neumann',
                  xaxis_title="x [m]",
                  yaxis_title="T [K]")


plot(fig, auto_open=True)
