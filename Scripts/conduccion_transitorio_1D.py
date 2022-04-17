# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 16:06:39 2022

@author: Rosendo
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot
from plotly.tools import make_subplots

# CONSTANTES ----------------------------------------------------------------#
L = 3        # Largo de barra [m]
S = 1        # Área transversal [m2]
Q = 0.01     # Fuente/(rho*Cp) [°C/s]
alpha = 0.01 # Difusividad [m2/s]  k/(rho*Cp)
n = 10       # 
t_max = 800    # Tiempo máximo en [s]

# Condiciones de Borde de Dirichlet (sobre la temperatura) ------------------#
phi_0 = 0
phi_L = 0

# Condicion inicial ---------------------------------------------------------#
phi_init = np.zeros(n)

#Theta 0: Forwawrd Euler | 1/2: Crank Nicolson | 1: Backward Euler 

def analytic_solution(x):
    c0 = phi_0
    c1 = phi_L/L + Q*L/(2*alpha) - c0/L
    return -Q*np.power(x,2)/(2*alpha) + c1*x + c0

def numerical_transitory_solution(n, phi_init, theta, t_max, deltaT, second_order_dirichlet = False):
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



# Add traces, one for each slider step
deltaT = 4.4
Fo_low = 2*alpha*deltaT/(L/n)**2
print('Fourier =', Fo_low)
theta = 0    # Paso temporal [s]
results_forward_lowFo = numerical_transitory_solution(n, phi_init, theta, t_max, deltaT, second_order_dirichlet = False)
theta = 0.5
results_crank_lowFo = numerical_transitory_solution(n, phi_init, theta, t_max, deltaT, second_order_dirichlet = False)
theta = 1
results_backward_lowFo = numerical_transitory_solution(n, phi_init, theta, t_max, deltaT, second_order_dirichlet = False)

deltaT = 4.7

theta = 0    # Paso temporal [s]
results_forward_highFo = numerical_transitory_solution(n, phi_init, theta, t_max, deltaT, second_order_dirichlet = False)
theta = 0.5
results_crank_highFo = numerical_transitory_solution(n, phi_init, theta, t_max, deltaT, second_order_dirichlet = False)
theta = 1
results_backward_highFo = numerical_transitory_solution(n, phi_init, theta, t_max, deltaT, second_order_dirichlet = False)

Fo_high = 2*alpha*deltaT/(L/n)**2
print('Fourier =', Fo_high)

fig = make_subplots(2, 3,
                    subplot_titles=['Forward Euler - Fo = ' + str(round(Fo_low, 4)),
                                    'Crank Nicolson - Fo = ' + str(round(Fo_low, 4)),
                                    'Backward Euler - Fo = ' + str(round(Fo_low, 4)),
                                    'Forward Euler - Fo = ' + str(round(Fo_high, 4)),
                                    'Crank Nicolson - Fo = ' + str(round(Fo_high, 4)),
                                    'Backward Euler - Fo = ' + str(round(Fo_high, 4))])

for time in results_forward_lowFo:
    if time == 'centroides':
        continue
    fig.add_scatter(
            row=1,
            col=1,
            visible=False,
            #line=dict(color="#00CED1", width=6),
            name="t = " + time + " [s]",
            x=np.concatenate(([0], results_forward_lowFo['centroides'], [L])),
            y=np.concatenate(([phi_0], results_forward_lowFo[time].values, [phi_L])))

for time in results_crank_lowFo:
    if time == 'centroides':
        continue
    fig.add_scatter(
            row=1,
            col=2,
            visible=False,
            #line=dict(color="#00CED1", width=6),
            name="t = " + time + " [s]",
            x=np.concatenate(([0], results_crank_lowFo['centroides'], [L])),
            y=np.concatenate(([phi_0], results_crank_lowFo[time].values, [phi_L])))

for time in results_backward_lowFo:
    if time == 'centroides':
        continue
    fig.add_scatter(
            row=1,
            col=3,
            visible=False,
            #line=dict(color="#00CED1", width=6),
            name="t = " + time + " [s]",
            x=np.concatenate(([0], results_backward_lowFo['centroides'], [L])),
            y=np.concatenate(([phi_0], results_backward_lowFo[time].values, [phi_L])))
    
for time in results_forward_highFo:
    if time == 'centroides':
        continue
    fig.add_scatter(
            row=2,
            col=1,
            visible=False,
            #line=dict(color="#00CED1", width=6),
            name="t = " + time + " [s]",
            x=np.concatenate(([0], results_forward_highFo['centroides'], [L])),
            y=np.concatenate(([phi_0], results_forward_highFo[time].values, [phi_L])))

for time in results_crank_highFo:
    if time == 'centroides':
        continue
    fig.add_scatter(
            row=2,
            col=2,
            visible=False,
            #line=dict(color="#00CED1", width=6),
            name="t = " + time + " [s]",
            x=np.concatenate(([0], results_crank_highFo['centroides'], [L])),
            y=np.concatenate(([phi_0], results_crank_highFo[time].values, [phi_L])))

for time in results_backward_highFo:
    if time == 'centroides':
        continue
    fig.add_scatter(
            row=2,
            col=3,
            visible=False,
            #line=dict(color="#00CED1", width=6),
            name="t = " + time + " [s]",
            x=np.concatenate(([0], results_backward_highFo['centroides'], [L])),
            y=np.concatenate(([phi_0], results_backward_highFo[time].values, [phi_L])))
    

len_lowFo = results_forward_lowFo.shape[1]
len_highFo = results_forward_highFo.shape[1]
# Make 10th trace visible
fig.data[0].visible = True
fig.data[len_lowFo-1].visible = True
fig.data[2*(len_lowFo-1)].visible = True
fig.data[3*(len_lowFo-1)].visible = True
fig.data[3*(len_lowFo-1)+len_highFo-1].visible = True
fig.data[3*(len_lowFo-1)+2*(len_highFo-1)].visible = True

# Create and add slider
steps = []
for i in range(len_highFo):
    step = dict(
        method="update",
        args=[{"visible": [False] * (len(fig.data))}]
    )
    step["args"][0]["visible"][i-1] = True
    step["args"][0]["visible"][i+len_lowFo-1-1] = True  
    step["args"][0]["visible"][i+2*(len_lowFo-1)-1] = True  
    step["args"][0]["visible"][i+3*(len_lowFo-1)-1] = True  
    step["args"][0]["visible"][i+3*(len_lowFo-1)+len_highFo-1-1] = True  
    step["args"][0]["visible"][i+3*(len_lowFo-1)+2*(len_highFo-1)-1] = True  
    steps.append(step)

sliders = [dict(
    active=0,
    currentvalue={"prefix": "Step: "},
    pad={"t": 50},
    steps=steps
)]

fig.update_layout(
    sliders=sliders,
    yaxis1_range=[-0.5, 1.5],
    yaxis2_range=[-0.5, 1.5],
    yaxis3_range=[-0.5, 1.5],
    yaxis4_range=[-0.5, 1.5],
    yaxis5_range=[-0.5, 1.5],
    yaxis6_range=[-0.5, 1.5]
)

plot(fig, auto_open=True)
