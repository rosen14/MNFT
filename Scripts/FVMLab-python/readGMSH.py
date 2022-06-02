# -*- coding: utf-8 -*-
"""
Created on Tue May 31 20:18:05 2022

@author: Rosendo
"""
import numpy as np

def fskipl(f, n):
    for i in range(n):
        next(f)
    return
    
def readGMSH(filename):
    '''
     Reads a GMSH .msh file
     Returns xnod and icone files
     with node coordinates and element conectivities
     Mesh generation example: gmsh -2 barra1D.geo -format msh2
     (C) 2017 - Santiago Marquez Damian
    '''
    
    # Open file for reading
    with open(filename, 'r') as fid:
        
        fskipl(fid, 4)

        # Get number of nodes
        nnod = int(fid.readline())
       
        # Alocate memory for nodes
        xnod = np.zeros((nnod, 3))
        
        # Get nodes coordinates
        for i in range(1,nnod + 1):
            # Read line by line
            # Scan for values with templated input
            sLine = [float(i) for i in fid.readline().split(' ')[1:]]
            
            # Store them discarding point number
            xnod[i-1] = np.array(sLine, dtype = 'float') #¿Precisión?
    
        # Avoid header
        fskipl(fid, 2)

        nel = int(fid.readline())
    
        #  Get element conectivities
        #  Some elements are lines at boundaries
        #  and points, intended to be at the beginning
        #  Mixed elements not supported
        for i in range(nel):
            # Read line by line
            ss = fid.readline()
            
            # Split line by spaces
            split = ss.split(' ')
           
            # Read element type
            elType = float(split[1]);
    
            # Continue reading until new type of elements
            if (elType != 1 and elType != 15):
                break
        
        nel = nel - i - 1;
    
        # Allocate memory for conectivities
        # Type of array depends on type of element
        if (elType == 2):
            # 3-node triangle
            icone = np.zeros((nel, 3))
            nnel = 3
        elif (elType == 3):
            # 4-node quadrangle
            icone = np.zeros((nel, 4))
            nnel = 4

        # Take number of data in splited line
        ndata = len(split)
        
        # Process last read line
        # Assert data number
        # elm-number elm-type number-of-tags < tag > … node-number-list
        
        # Take number of tags
        ntags = int(split[2])
        
        if (ndata - 3 - ntags != nnel):
            print('ERROR reading GMSH file 1\n');
            return
        
        # Alocate memory for elements
        icone = np.zeros((nel + 1, nnel))
        
        # Take data from last read line
        j = 0
        for i in range((3 + ntags), ndata):
            icone[0, j] = float(split[i])
            j+=1
        
        # Read conectivities
        for i in range(nel):
            # Read line by line
            ss = fid.readline()
            
            # Split line by spaces
            split = ss.split(' ')
            
            # Take number of data in splited line
            ndata = len(split)
            
            if (ndata - 3 - ntags != nnel):
                print('ERROR reading GMSH file\n');
                return
            
            # Take data from last read line
            k = 0
            for j in range((3 + ntags), ndata):
                icone[i+1, k] = float(split[j])
                k+=1
    return xnod, icone