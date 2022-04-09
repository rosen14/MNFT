# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 20:56:33 2022

@author: Rosendo
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot

L = 3    # Largo de barra [m]
n = 10   # Discretización del dominio
h = L/n  # Tamaño de discretizacion (uniforme)
S = 1    # Área transversal [m2]
Q = 0.01 # Fuente [W/m2]
k = 0.01 # Conductividad [W/mK]

#Condiciones de Borde de Dirichlet (sobre la temperatura)
phi_0 = 0
phi_L = 0

def analytic_solution(x, phi_0, phi_L, L, Q, k):
    c0 = phi_0
    c1 = phi_L/L + Q*L/(2*k) - c0/L
    return -Q*np.power(x,2)/(2*k) + c1*x + c0

# Coordenada de caras
x_caras = np.linspace(0, L, n + 1)

#Coordenada de centroide de las celdas
x_centroides = np.array([0.5*(x_caras[i+1] - x_caras[i]) + x_caras[i] \
              for i in range(len(x_caras) - 1)])

b_matrix = np.zeros((n,))
A_matrix = np.zeros((n,n))

for i in range(1, n-1):
    A_matrix[i,i-1:i+2] = np.array([1, -2, 1])
    b_matrix[i] = -Q*h**2/k
    
#Por condición de borde
A_matrix[0, 0:2] = np.array([-3, 1])
A_matrix[n-1, n-2:n] = np.array([1, -3])

b_matrix[0] = -Q*h**2/k - 2*phi_0
b_matrix[n-1] = -Q*h**2/k - 2*phi_L

phi_matrix = np.linalg.solve(A_matrix, b_matrix)

reference_sol = [0.225000000000000,
                 0.585000000000001,
                 0.855000000000001,
                 1.035000000000002,
                 1.125000000000002,
                 1.125000000000002,
                 1.035000000000002,
                 0.855000000000001,
                 0.585000000000001,
                 0.225000000000000]

if n == 10:
    diff_df = pd.DataFrame({'x': x_centroides,
                            'Con referencia': phi_matrix - reference_sol,
                            'Con analitica': phi_matrix - analytic_solution(
                               x_centroides, phi_0, phi_L, L, Q, k)
                            })

##PLOT###########
fig = go.Figure()

# Add traces
fig.add_trace(
    go.Scatter(
        x = np.concatenate(([0], x_centroides, [L])),
        y = np.concatenate(([phi_0], phi_matrix, [phi_L])),
        mode='markers', name = 'VF',
        marker=dict(
        size=8))
    )

analytic_domain = np.linspace(0, L, 100)
fig.add_trace(
    go.Scatter(
        x = analytic_domain,
        y = analytic_solution(analytic_domain, phi_0, phi_L, L, Q, k),
        mode='lines+markers', name = 'Analytic'))

fig.update_layout(title='Difusión estacionario 1D',
                  xaxis_title="x [m]",
                  yaxis_title="T [K]")


plot(fig, auto_open=True)
#################