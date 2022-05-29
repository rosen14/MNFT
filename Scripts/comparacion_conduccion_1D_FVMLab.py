# -*- coding: utf-8 -*-
"""
Created on Sun May 29 19:37:38 2022

@author: Rosendo
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot

# CONSTANTES ----------------------------------------------------------------#
L = 1        # Largo de barra [m]
S = 1        # Área transversal [m2]
Q = 0.01  # Fuente/(rho*Cp) [°C/s]
alpha = 0.5 # Difusividad [m2/s]  k/(rho*Cp)


# Condiciones de Borde de Dirichlet (sobre la temperatura) ------------------#
phi_0 = 0
phi_L = 0

def analytic_solution(x):
    c0 = phi_0
    c1 = phi_L/L + Q*L/(2*alpha) - c0/L
    return -Q*np.power(x,2)/(2*alpha) + c1*x + c0


def numerical_solution(n, second_order_dirichlet = False):
    h = L/n     # Tamaño de discretizacion (uniforme)
    
    x_caras = np.linspace(0, L, n + 1)  # Coordenada de caras
    
    # Coordenada de centroide de las celdas
    x_centroides = np.array([0.5*(x_caras[i+1] - x_caras[i]) + x_caras[i] \
                  for i in range(len(x_caras) - 1)])
    
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
        
    phi_numerical = np.linalg.solve(K_matrix, b_matrix)

    return x_centroides, phi_numerical


# COMPARACIÓN CON SOLUCION REFERENCIA ---------------------------------------#
xc_n10, phi_n10 = numerical_solution(n = 10, second_order_dirichlet=False)

FVMLab_sol = [5.0000e-04,
   1.3000e-03,
   1.9000e-03,
   2.3000e-03,
   2.5000e-03,
   2.5000e-03,
   2.3000e-03,
   1.9000e-03,
   1.3000e-03,
   5.0000e-04]
diff_df = pd.DataFrame({'x': xc_n10,
                        'Con referencia': phi_n10 - FVMLab_sol,
                        })
#----------------------------------------------------------------------------#

# PLOT ----------------------------------------------------------------------#
fig = go.Figure()

x_domain = np.concatenate(([0], xc_n10, [L]))

# Add traces
fig.add_trace(
    go.Scatter(
        x = x_domain,
        y = np.concatenate(([phi_0], phi_n10, [phi_L])),
        mode='lines+markers', name = 'Numerical',
        marker=dict(
        size=8))
    )
fig.add_trace(
    go.Scatter(
        x = x_domain,
        y = np.concatenate(([phi_0], FVMLab_sol, [phi_L])),
        mode='lines+markers', name = 'FVMLab_sol'))

fig.update_layout(title='Difusión estacionario 1D',
                  xaxis_title="x [m]",
                  yaxis_title="T [K]")


plot(fig, auto_open=True)