# -*- coding: utf-8 -*-
"""
Created on Sun May 15 22:12:32 2022

@author: Rosendo
"""


import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot

# CONSTANTES ----------------------------------------------------------------#
L = 1        # Largo de barra [m]
S = 1        # Área transversal [m2]
Q = 0     # Fuente/(rho*Cp) [°C/s]
alpha = 0.05 # Difusividad [m2/s]  k/(rho*Cp)
V = 1
# Condiciones de Borde de Dirichlet (sobre la temperatura) ------------------#
phi_0 = 0
phi_L = 1

def analytic_solution(x):
    PeL = 0.5*V*L/alpha
    return (1 - np.exp(2*PeL*x/L))/(1 - np.exp(2*PeL))


def numerical_solution(n, upwind = True, second_order_dirichlet = True):
    h = L/n     # Tamaño de discretizacion (uniforme)
    
    x_caras = np.linspace(0, L, n + 1)  # Coordenada de caras
    
    # Coordenada de centroide de las celdas
    x_centroides = np.array([0.5*(x_caras[i+1] - x_caras[i]) + x_caras[i] \
                  for i in range(len(x_caras) - 1)])
    
    b_matrix = np.zeros((n,))  #Matriz término fuente y C.B.
    K_matrix = np.zeros((n,n)) #Matriz de difusión
    
    ba_matrix = np.zeros((n,))
    A_matrix = np.zeros((n,n)) #Matriz de advección
    
    if second_order_dirichlet:
        K_matrix[0, 0:2] = -alpha/h*np.array([-4, 4/3])
        K_matrix[n-1, n-2:n] = -alpha/h*np.array([4/3, -4])
        
        b_matrix[0] = Q*h + 8/3*alpha/h*phi_0
        b_matrix[n-1] = Q*h + 8/3*alpha/h*phi_L
    else:
        # Por condición de borde
        K_matrix[0, 0:2] = -alpha/h*np.array([-3, 1])
        K_matrix[n-1, n-2:n] = -alpha/h*np.array([1, -3])
        
        b_matrix[0] = Q*h + 2*alpha/h*phi_0
        b_matrix[n-1] = Q*h + 2*alpha/h*phi_L
            
        #Advección
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


    for i in range(1, n-1):
        K_matrix[i,i-1:i+2] = -alpha/h*np.array([1, -2, 1])
        b_matrix[i] = Q*h
        
    
    phi_numerical = np.linalg.solve(K_matrix + A_matrix, b_matrix + ba_matrix)

    return x_centroides, phi_numerical


#----------------------------- Tarea 19/05/2022 -----------------------------#
## a) PeL = 10 ; N = 40 celdas
#### i) Calcular nro de Pe de malla: 
n = 40
h = L/n 
Pex = 0.5*V*h/alpha
print('Pe_x = ', Pex)
print('Pe_L = ', 0.5*V*L/alpha)

#### ii) Solución con central differences
x_centroides, phi_numerical_cd = numerical_solution(n, upwind = False, second_order_dirichlet = True)
#### iii) Solución con upwind
x_centroides, phi_numerical_uw = numerical_solution(n, upwind = True, second_order_dirichlet = True)

#### iv) Comparar gráficamente CD, Upwind y la solución analítica calculada en
####     200 puntos
x_domain = np.linspace(0, L, 200)
phi_analytic = analytic_solution(x_domain)

fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x = np.concatenate(([0], x_centroides, [L])),
        y = np.concatenate(([phi_0], phi_numerical_cd, [phi_L])),
        mode='lines+markers', name = 'Central Differences',
        marker=dict(
        size=8))
    )
fig.add_trace(
    go.Scatter(
        x = np.concatenate(([0], x_centroides, [L])),
        y = np.concatenate(([phi_0], phi_numerical_uw, [phi_L])),
        mode='lines+markers', name = 'Upwind',
        marker=dict(
        size=8))
    )
fig.add_trace(
    go.Scatter(
        x = x_domain,
        y = phi_analytic,
        mode='lines', name = 'Analytic'))

plot(fig, auto_open=True)

#-------------

## b) PeL = 5 ; N = 8 celdas
#### i) Calcular nro de Pe de malla: 
n = 8
h = L/n 
Pex = 0.5*V*h/alpha
print('Pe_x = ', Pex)
print('Pe_L = ', 0.5*V*L/alpha)

#### ii) Solución con central differences
x_centroides, phi_numerical_cd = numerical_solution(n, upwind = False, second_order_dirichlet = False)
#### iii) Solución con upwind
x_centroides, phi_numerical_uw = numerical_solution(n, upwind = True, second_order_dirichlet = False)

#### iv) Comparar gráficamente CD, Upwind y la solución analítica calculada en
####     200 puntos
x_domain = np.linspace(0, L, 200)
phi_analytic = analytic_solution(x_domain)

fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x = np.concatenate(([0], x_centroides, [L])),
        y = np.concatenate(([phi_0], phi_numerical_cd, [phi_L])),
        mode='lines+markers', name = 'Central Differences',
        marker=dict(
        size=8))
    )
fig.add_trace(
    go.Scatter(
        x = np.concatenate(([0], x_centroides, [L])),
        y = np.concatenate(([phi_0], phi_numerical_uw, [phi_L])),
        mode='lines+markers', name = 'Upwind',
        marker=dict(
        size=8))
    )
fig.add_trace(
    go.Scatter(
        x = x_domain,
        y = phi_analytic,
        mode='lines', name = 'Analytic'))

plot(fig, auto_open=True)

## c) Analisis de convergencia del error con PeL = 5 ; N = 10 celdas de base,
##    i) Con upwind, resultados en log-log
##    ii) Con Central differences, resultados en log-log

alpha = 0.1 #Para PeL = 5
print('Pe_x = ', 0.5*V*h/alpha)
print('Pe_L = ', 0.5*V*L/alpha)

err_N = [10*2**(i-1) for i in range(1, 8)]


error_df = pd.DataFrame()
for N in err_N:
    xc, phi_numerical_cd = numerical_solution(N, upwind = False, second_order_dirichlet = False)
    xc, phi_numerical_uw = numerical_solution(N, upwind = True, second_order_dirichlet = False)
    phi_analytic = analytic_solution(xc)
    
    norma1_cd = 1/N * np.sum(np.abs(phi_analytic - phi_numerical_cd))
    norma1_uw = 1/N * np.sum(np.abs(phi_analytic - phi_numerical_uw))
    
    norma2_cd = np.sqrt(1/N * np.sum(np.abs(phi_analytic - phi_numerical_cd)**2))
    norma2_uw = np.sqrt(1/N * np.sum(np.abs(phi_analytic - phi_numerical_uw)**2))
    
    error_df = error_df.append({
                    'discretizacion': N,
                    'norma 1 - CD': norma1_cd,
                    'norma 2 - CD': norma2_cd,
                    'norma 1 - UW': norma1_uw,
                    'norma 2 - UW': norma2_uw,
                    },
                    ignore_index=True)


error_df['tasa error n1 - CD'] = pd.Series(np.log2(error_df['norma 1 - CD'][:-1].values/ \
                                      error_df['norma 1 - CD'][1:].values))

error_df['tasa error n1 - UW'] = pd.Series(np.log2(error_df['norma 1 - UW'][:-1].values/ \
                                      error_df['norma 1 - UW'][1:].values))
    
error_df['tasa error n2 - CD'] = pd.Series(np.log2(error_df['norma 2 - CD'][:-1].values/ \
                                      error_df['norma 2 - CD'][1:].values))

error_df['tasa error n2 - UW'] = pd.Series(np.log2(error_df['norma 2 - UW'][:-1].values/ \
                                      error_df['norma 2 - UW'][1:].values))

fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x = err_N,
        y = error_df['norma 1 - CD'],
        mode='lines+markers', name = 'Norma 1 - CD',
        marker=dict(
        size=8))
    )
fig.add_trace(
    go.Scatter(
        x = err_N,
        y = error_df['norma 2 - CD'],
        mode='lines+markers', name = 'Norma 2 - CD',
        marker=dict(
        size=8))
    )
fig.add_trace(
    go.Scatter(
        x = err_N,
        y = error_df['norma 1 - UW'],
        mode='lines+markers', name = 'Norma 1 - UW',
        marker=dict(
        size=8))
    )
fig.add_trace(
    go.Scatter(
        x = err_N,
        y = error_df['norma 2 - UW'],
        mode='lines+markers', name = 'Norma 2 - UW',
        marker=dict(
        size=8))
    )
fig.update_xaxes(type="log")
fig.update_yaxes(type="log")
plot(fig, auto_open=True)


