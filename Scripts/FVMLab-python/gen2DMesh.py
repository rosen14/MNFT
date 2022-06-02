# -*- coding: utf-8 -*-
"""
Created on Tue May 31 20:14:27 2022

@author: Rosendo
"""

class mesh:
    def __init__(self):
        self.xnod = None
        self.icone = None
        self.C = None
        
        self.mesh.faces = []
        self.mesh.owner = []
        self.mesh.neighbour = []
        self.mesh.Sf = []
        self.mesh.Cf = []
        self.mesh.V = []
        self.mesh.iface = []
    
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

    mesh = mesh()
    # GMSH reading
    mesh.xnod, mesh.icone = readGMSH(filename)
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
            f = [mesh.icone[i,j], mesh.icone[i, (j % mesh.nen) + 1]]
            
            # Check if this face already exist
            mesh.iface = []
            if mesh.faces:
                mesh.iface = np.where(mesh.faces[:, 0] == f[1] & mesh.faces[:, 1] == f[1])
            end
            if (isempty(mesh.iface))
                mesh.faces = [mesh.faces; f];
                mesh.owner = [mesh.owner; i];
                mesh.neighbour = [mesh.neighbour; 0];
                vf = mesh.xnod(f(2),:) - mesh.xnod(f(1),:);
                cf = (mesh.xnod(f(1),:) + mesh.xnod(f(2),:))/2.0;
                mesh.Sf = [mesh.Sf; vf(2),-vf(1)];
                mesh.Cf = [mesh.Cf; cf];   
                %quiver(cf(1),cf(2),mesh.Sf(end,1),mesh.Sf(end,2),0.25,'b');
                v = v + cf*mesh.Sf(end, :)';
            else
                mesh.neighbour(mesh.iface) = i;
                v = v + mesh.Cf(mesh.iface, :)*(-1)*mesh.Sf(mesh.iface, :)';
            end
            % Volumes from Gauss theorem
        end
        mesh.V = [mesh.V; v/2];
    end
    
        
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