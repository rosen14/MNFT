# -*- coding: utf-8 -*-
"""
Created on Sun May  8 22:14:00 2022

@author: Rosendo
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import plot
import os

# CONSTANTES ----------------------------------------------------------------#
L = 3        # Largo de barra [m]
S = 1        # Área transversal [m2]
Q = 0.01     # Fuente/(rho*Cp) [°C/s]
alpha = 0.01 # Difusividad [m2/s]  k/(rho*Cp)

# Condiciones de Borde de Dirichlet (sobre la temperatura) ------------------#
phi_0 = 0
phi_L = 0

def analytic_solution(x):
    c0 = phi_0
    c1 = phi_L/L + Q*L/(2*alpha) - c0/L
    return -Q*np.power(x,2)/(2*alpha) + c1*x + c0



error_df = pd.DataFrame()
folder_path = 'Convergencia_Malla_DiffEst1D_OpemFOAM/1/'
for filename in os.listdir(folder_path):
    open_foam_sol = pd.read_csv(os.path.join(folder_path, filename),
                    sep= '\t', names=['xc', 'T'])
    
    N = len(open_foam_sol)
    
    phi_analytic = analytic_solution(open_foam_sol['xc'])
    
    norma1 = 1/N * np.sum(np.abs(phi_analytic - open_foam_sol['T']))
    norma2 = np.sqrt(1/N * np.sum(np.abs(phi_analytic - open_foam_sol['T'])**2))
    
    error_df = error_df.append({
                    'discretizacion': N,
                    'norma 1': norma1,
                    'norma 2': norma2
                    },
                    ignore_index=True)
    
error_df = error_df.sort_values(by=['discretizacion']).reset_index(drop=True)

error_df['tasa error n1'] = pd.Series(np.log2(error_df['norma 1'][:-1].values/ \
                                      error_df['norma 1'][1:].values))

error_df['tasa error n2'] = pd.Series(np.log2(error_df['norma 2'][:-1].values/ \
                                      error_df['norma 2'][1:].values))
    