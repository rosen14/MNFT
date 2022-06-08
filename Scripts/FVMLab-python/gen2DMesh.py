# -*- coding: utf-8 -*-
"""
Created on Tue May 31 20:14:27 2022

@author: Rosendo
"""

import numpy as np
import readGMSH
class Mesh:
    def __init__(self):
        self.xnod = None
        self.icone = None
        self.C = None
        
        self.faces = []
        self.owner = []
        self.neighbour = []
        self.Sf = []
        self.Cf = []
        self.V = []
        self.iface = []

    
def gen2DMesh(filename):
    '''
     Generates 2D mesh information reading
    from GMSH mesh format

    (C) 2016 - 2017 - Juan Marcelo Gimenez
    (C) 2017 - Santiago Marquez Damian
                - Formatting
                - Native GMSH data support
                - Element to face map added (mesh.ele)
    '''

    mesh = Mesh()
    # GMSH reading
    mesh.xnod, mesh.icone = readGMSH.readGMSH(filename)
    # Override third dimension in points
    mesh.xnod = mesh.xnod[:, 0:2]

    mesh.nelem, mesh.nen = mesh.icone.shape
    mesh.nnod, mesh.ndim = mesh.xnod.shape

    # Centroids
    mesh.C = np.zeros((mesh.nelem,2))
    for i in range(mesh.nen):
        mesh.C[:,0] = mesh.C[:,0] + np.take(mesh.xnod[:,0], [mesh.icone[:,i] - 1])
        mesh.C[:,1] = mesh.C[:,1] + np.take(mesh.xnod[:,1], [mesh.icone[:,i] - 1])
        
    mesh.C = mesh.C/mesh.nen


    delta = []
    k = []
    
    for i in range(mesh.nelem):
        v = 0
        for j in range(mesh.nen):
            f = [int(mesh.icone[i,j]), int(mesh.icone[i, (j % mesh.nen)])]
            
            # Check if this face already exist
            mesh.iface = -1
            if len(mesh.faces):
                mesh.iface = np.where(np.array(mesh.faces)[:, 0] == f[1] and np.array(mesh.faces)[:, 1] == f[1])[0][0]
                
            if mesh.iface == -1:
                mesh.faces.append(f)
                mesh.owner.append(i)
                mesh.neighbour.append(0)
                vf = mesh.xnod[f[1],:] - mesh.xnod[f[0],:]
                cf = (mesh.xnod[f[0],:] + mesh.xnod[f[1],:])/2.0
                mesh.Sf.append([vf[1],-vf[0]])
                mesh.Cf.append(list(cf))   
                
                v = v + np.dot(cf, mesh.Sf[-1])
            else:
                print('a')
                mesh.neighbour[mesh.iface] = i;
                v = v - np.dot(np.take(mesh.Cf, mesh.iface), np.take(mesh.Sf, mesh.iface))

        mesh.V.append(v/2)
    
        
    % To store element to face map
    % -> Could this be done in previous loop?
    % First column is pointer
    mesh.ele = zeros(mesh.nelem, mesh.nen + 1);
    mesh.ele(:, 1) = 2;
    
    % Walk along owner and neighbour adding them to cells
    for i = 1:size(mesh.owner)
        % Take owner and neighbour for present face
        owner = mesh.owner(i);
        neighbour = mesh.neighbour(i);
        
        % Take pointer and store face number
        j = mesh.ele(owner, 1);
        mesh.ele(owner, j) = i;
        
        % Increment pointer
        mesh.ele(owner, 1) = mesh.ele(owner, 1) + 1;
        
        % If neighbour isn't zero store it
        if (neighbour != 0)
            % Take pointer and store face number
            j = mesh.ele(neighbour, 1);
            mesh.ele(neighbour, j) = i;
            
            % Increment pointer
            mesh.ele(neighbour, 1) = mesh.ele(neighbour, 1) + 1;
        end
    end
    
    % Drop first column
    mesh.ele = mesh.ele(:, 2:end);
    
    % Alias
    mesh.ncells = mesh.nelem;
    
    % Number of faces
    mesh.nfaces = size(mesh.faces, 1);
return mesh