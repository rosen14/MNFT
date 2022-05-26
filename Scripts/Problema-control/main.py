# -*- coding: utf-8 -*-
"""
Created on Sat May 21 17:35:35 2022

@author: Rosendo
"""


#import plotly.graph_objects as go
#from plotly.offline import plot
import numpy as np
import scipy.integrate as integrate
import pandas as pd
from scipy.sparse import lil_matrix, identity
from scipy.sparse.linalg import cgs, spilu, LinearOperator
from math import log10, floor


def ss_sol(nu, h, phiInf, x, L):
    b = np.sqrt(h/(L*nu))
    sol = phiInf*(1 - np.cosh(b*x)) + \
          phiInf*(np.cosh(b*L) - 1)*np.sinh(b*x)/np.sinh(b*L)
    return sol

def lambda_n_2(nu, h, L, n):
    return (n*np.pi/L)**2*nu + h/L

def integrando(x, nu, h, phiInf, L, n):
    return ss_sol(nu, h, phiInf, x, L)*np.sin(n*np.pi*x/L)

def cn(nu, h, phiInf, L, n):
    cn = -2*integrate.quad(integrando, 0, L, args=(nu, h, phiInf, L, n))[0]/L
    return cn


def TranReacDiff(nu, h, phiInf, x, L, t, epsilon, nmax):
    phi = ss_sol(nu, h, phiInf, x, L)
    update_phi = phi.copy()
    N = len(x)
    for n in np.arange(1, nmax+1):
        update_phi += cn(nu, h, phiInf, L, n)* \
               np.sin(n*np.pi*x/L)* \
               np.exp(-lambda_n_2(nu, h, L, n)*t)
               
        norma1 = 1/N * np.sum(np.abs(update_phi - phi))
        phi = update_phi.copy()
        if norma1 <= epsilon and n%2:
            # print('epsilon')
            # print(n)
            break
    return phi

def numerical_transitory_solution(n, phi_init, theta, t_max, deltaT, 
                                  nu, h, phiInf, phi_0, phi_L, L,
                                  second_order_dirichlet = False):
    l = L/n     # Tamaño de discretizacion (uniforme)

    x_caras = np.linspace(0, L, n + 1)  # Coordenada de caras
    
    # Coordenada de centroide de las celdas
    x_centroides = np.array([0.5*(x_caras[i+1] - x_caras[i]) + x_caras[i] \
                  for i in range(len(x_caras) - 1)])
    
    df_transitory = pd.DataFrame({'centroides': x_centroides})
    
    b_matrix = np.zeros((n,))  #Matriz término fuente y C.B.
    K_matrix = lil_matrix((n,n))
    
    # Por condición de borde
    if second_order_dirichlet:
        K_matrix[0, 0:2] = np.array([4*nu/l + h*l/L, -4*nu/(3*l)])
        K_matrix[n-1, n-2:n] = np.array([-4*nu/(3*l), 4*nu/l + h*l/L])
        
        b_matrix[0] = h*phiInf*l/L + 8/3*nu/l*phi_0
        b_matrix[n-1] = h*phiInf*l/L + 8/3*nu/l*phi_L
    else:
        K_matrix[0, 0:2] = np.array([3*nu/l + h*l/L, -nu/l])
        K_matrix[n-1, n-2:n] = np.array([-nu/l, 3*nu/l + h*l/L])
        
        b_matrix[0] = h*phiInf*l/L + 2*nu/l*phi_0
        b_matrix[n-1] = h*phiInf*l/L + 2*nu/l*phi_L

    for i in range(1, n-1):
        K_matrix[i,i-1:i+2] = np.array([-nu/l, 2*nu/l + h*l/L, -nu/l])
        b_matrix[i] = h*phiInf*l/L
        
    K_matrix = K_matrix.tocsc()
    I_matrix = identity(n).tocsc()
    M_matrix = l/deltaT*I_matrix + K_matrix*theta
    
    sA_iLU = spilu(M_matrix, drop_tol = 1e-3)
    LU = LinearOperator((n,n), sA_iLU.solve)
    
    phi_k = phi_init
    t = 0
    steps = int(t_max/deltaT)
    for i in range(steps):
        
        R_matrix = b_matrix + phi_k/deltaT*l - (1 - theta)*K_matrix.dot(phi_k)
    
        phi_k1, flag = cgs(M_matrix, R_matrix, x0 = phi_k, tol = 1e-9, M = LU)
              
        phi_k = phi_k1
        if flag:
            print('Sin convergencia, deteniendose...')
            return 0
            
    df_transitory[str(t_max)] = phi_k1
    
    return df_transitory

def get_convergence_order(ndiv):
    phi_0 = 0
    phi_L = 0
    nu = 1 #Difusividad
    L = 3
    phiInf = 10
    epsilon = 1e-12
    nmax = 100
    
    x_caras = np.linspace(0, L, ndiv + 1)  # Coordenada de caras
    x_centroides = np.array([0.5*(x_caras[i+1] - x_caras[i]) + x_caras[i] \
                  for i in range(len(x_caras) - 1)])
        
    phi_init = np.zeros(ndiv)
    
    dTBase = 0.01
    discr_temp = [dTBase/2**(i-1) for i in range(1, 6)] 
    t_max = 0.1
    
    # Solucion de referencia ####
    h = 10
    
    
    ref_solution = TranReacDiff(nu, h, phiInf, x_centroides, L, t_max, epsilon, nmax)
    ###############################
    rangeIndex = pd.RangeIndex(start=1, stop=len(discr_temp), step=1) #Para indexar bien la tasa de error
    error_df = pd.DataFrame({'Paso de tiempo': [],
                             'Fo': [],
                             'error L1 BE': [],
                             'orden de convergencia BE': [],
                             'error L1 CN': [],
                             'orden de convergencia CN': []
                             })
    for deltaT in discr_temp:
        phi_numerical_crank = numerical_transitory_solution(ndiv, phi_init, 0.5,
                                                            t_max, deltaT, nu, h,
                                                            phiInf, phi_0, phi_L,
                                                            L, second_order_dirichlet = True)
        phi_numerical_backward = numerical_transitory_solution(ndiv, phi_init, 1,
                                                            t_max, deltaT, nu, h,
                                                            phiInf, phi_0, phi_L,
                                                            L, second_order_dirichlet = True)
        
        norma1_crank = 1/ndiv * np.sum(np.abs(ref_solution - phi_numerical_crank[str(t_max)]))
        norma1_backward = 1/ndiv * np.sum(np.abs(ref_solution - phi_numerical_backward[str(t_max)]))
        
        error_df = error_df.append({
                        'Paso de tiempo':deltaT,
                        'Fo': 2*nu*deltaT/(L/ndiv)**2,
                        'error L1 CN': norma1_crank,
                        'error L1 BE': norma1_backward
                        },
                        ignore_index=True)
        
    
    error_df['orden de convergencia CN'] = pd.Series(np.log2(error_df['error L1 CN'][:-1].values/ \
                                                          error_df['error L1 CN'][1:].values), index = rangeIndex)
    error_df['orden de convergencia BE'] = pd.Series(np.log2(error_df['error L1 BE'][:-1].values/ \
                                                          error_df['error L1 BE'][1:].values), index = rangeIndex)

    return error_df.fillna('-')

if __name__ == "__main__":
    error_df_ndiv_10000 = get_convergence_order(10000)
    error_df_ndiv_1000 = get_convergence_order(1000)
    
    pd.set_option('display.float_format', '{:.3E}'.format)
    
    print('10.000 elementos: \n ', error_df_ndiv_10000)
    print('\n 1.000 elementos: \n ', error_df_ndiv_1000)
    error_df_ndiv_1000


#------------------------------------------
'''
    
#Item d)
phi_0 = 0
phi_L = 0
nu = 1 #Difusividad
L = 3
phiInf = 10
epsilon = 1e-12
nmax = 100
t = 0.1

ndiv = 10000
x_caras = np.linspace(0, L, ndiv + 1)  # Coordenada de caras
x_centroides = np.array([0.5*(x_caras[i+1] - x_caras[i]) + x_caras[i] \
              for i in range(len(x_caras) - 1)])
    


sol_h10 = TranReacDiff(nu, 10, phiInf, x_centroides, L, t, epsilon, nmax)
sol_h100 = TranReacDiff(nu, 100, phiInf, x_centroides, L, t, epsilon, nmax)
sol_ss_h10 = ss_sol(nu, 10, phiInf, x_centroides, L)
sol_ss_h100 = ss_sol(nu, 100, phiInf, x_centroides, L)


###############################################################
#PLOTEO##

# Create figure
fig = go.Figure()

fig.add_trace(
        go.Scatter(
            name = 'h = 10 | t = 0.1',
            x=x_centroides,
            y=sol_h10))

fig.add_trace(
        go.Scatter(
            name = 'h = 100 | t = 0.1',
            x=x_centroides,
            y=sol_h100))

fig.add_trace(
        go.Scatter(
            name = 'h = 10 | t = inf',
            x=x_centroides,
            y=sol_ss_h10))

fig.add_trace(
        go.Scatter(
            name = 'h = 100 | t = inf',
            x=x_centroides,
            y=sol_ss_h100))

fig.update_layout(xaxis_title="x",
                  yaxis_title="phi",
                  autosize=False,
                  width=500,
                  height=500)

plot(fig, auto_open=True)

#-----------------

#Verificación mediante comparación con solución de referencia

n = 10000
phi_init = np.zeros(n)
theta = 0.5
t_max = 0.05
deltaT = 0.001
h = 10

cn_own = numerical_transitory_solution(n, phi_init, theta, t_max, deltaT, 
                                  nu, h, phiInf, phi_0, phi_L, L,
                                  second_order_dirichlet = True)

ref = np.loadtxt('phiRefNum.dat')

print('error contra referencia: ', 1/n * np.sum(np.abs(cn_own['0.05'] - ref)))

#---------------

dTBase = 0.01
discr_temp = [dTBase/2**(i-1) for i in range(1, 6)] 
nTdisc = len(discr_temp)
fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x = discr_temp,
        y = error_df_ndiv_10000['error L1 CN'],
        mode='lines+markers', name = 'error L1 CN - 1E4 elem.',
        marker=dict(
        size=8))
    )
fig.add_trace(
    go.Scatter(
        x = discr_temp,
        y = error_df_ndiv_10000['error L1 BE'],
        mode='lines+markers', name = 'error L1 BE - 1E4 elem.',
        marker=dict(
        size=8))
    )
fig.add_trace(
    go.Scatter(
        x = discr_temp,
        y = error_df_ndiv_1000['error L1 CN'],
        mode='lines+markers', name = 'error L1 CN - 1E3 elem.',
        marker=dict(
        size=8))
    )
fig.add_trace(
    go.Scatter(
        x = discr_temp,
        y = error_df_ndiv_1000['error L1 BE'],
        mode='lines+markers', name = 'error L1 BE - 1E3 elem.',
        marker=dict(
        size=8))
    )
fig.add_trace(
    go.Scatter(
        x = np.array([discr_temp[-1], discr_temp[0]]),
        y = np.array([error_df_ndiv_10000['error L1 BE'][nTdisc-1],
                      error_df_ndiv_10000['error L1 BE'][nTdisc-1]*np.power(2, 4)]),
        mode='lines', name = 'Orden 1 - Referencia',
        line=dict(color='black', dash='dot'))
    )

fig.add_trace(
    go.Scatter(
        x = np.array([discr_temp[-1], discr_temp[0]]),
        y = np.array([error_df_ndiv_10000['error L1 CN'][nTdisc-1],
                      error_df_ndiv_10000['error L1 CN'][nTdisc-1]*np.power(4, 4)]),
        mode='lines', name = 'Orden 2 - Referencia',
        line=dict(color='maroon', dash='dot'))
    )

fig.update_xaxes(type="log")
fig.update_yaxes(type="log")
fig.update_layout(xaxis_title="Discretización temporal",
                  yaxis_title="Error Norma 1",
                  autosize=False,
                  width=800,
                  height=500)

plot(fig, auto_open=True)
'''
