function [xnod, icone] = readGMSH(filename)
    % Reads a GMSH .msh file
    % Returns xnod and icone files
    % with node coordinates and element conectivities
    % Mesh generation example: gmsh -2 barra1D.geo -format msh2
    % (C) 2017 - Santiago Marquez Damian
    
    % Open file for reading
    fid = fopen(filename);    
    
    % Avoid header
    ss = fskipl(fid, 4);

    % Get number of nodes
    ss = fgetl(fid);
    nnod = str2num(ss);
    
    % Alocate memory for nodes
    xnod = zeros(nnod, 3);
    
    % Get nodes coordinates
    for i = 1:nnod
        % Read line by line
        ss = fgetl(fid);
        
        % Scan for values with templated input
        sLine = sscanf(ss, "%d %f %f %f");
        
        % Store them discarding point number
        xnod(i, 1:3) = sLine(2:4);
    end
    
    % Avoid header
    ss = fskipl(fid, 2);

    % Get number of elements
    ss = fgetl(fid);
    nel = str2num(ss);
    
    % Get element conectivities
    % Some elements are lines at boundaries
    % and points, intended to be at the beginning
    % Mixed elements not supported
    for i = 1:nel
        % Read line by line
        ss = fgetl(fid);
        
        % Split line by spaces
        split = strsplit(ss, ' ');
       
        % Read element type
        elType = str2num(split{2});

        % Continue reading until new type of elements
        if (elType != 1 && elType != 15)
            break;
        end
    end
    
    nel = nel - i;
    
    % Allocate memory for conectivities
    % Type of array depends on type of element
    if (elType == 2)
        % 3-node triangle
        icone = zeros(nel, 3);
        nnel = 3;
    elseif (elType == 3) 
        % 4-node quadrangle
        icone = zeros(nel, 4);
        nnel = 4;
    end
    
    % Take number of data in splited line
    ndata = size(split, 2);
    
    % Process last read line
    % Assert data number
    % elm-number elm-type number-of-tags < tag > … node-number-list
    
    % Take number of tags
    ntags = str2num(split{3});
    
    if (ndata - 3 - ntags != nnel)
        printf('ERROR reading GMSH file 1\n');
        return;
    end
    
    % Alocate memory for elements
    icone = zeros(nel + 1, nnel);
    
    % Take data from last read line
    j = 0;
    for i = (3 + ntags + 1):ndata
        j++;
        icone(1, j) = str2num(split{i});
    end
    
    % Read conectivities
    for i = 1:nel
        % Read line by line
        ss = fgetl(fid);
        
        % Split line by spaces
        split = strsplit(ss, ' ');
        
        % Take number of data in splited line
        ndata = size(split, 2);
        
        if (ndata - 3 - ntags != nnel)
            printf('ERROR reading GMSH file\n');
            return;
        end
        
        % Take data from last read line
        k = 0;
        for j = (3 + ntags + 1):ndata
            k++;
            icone(i + 1, k) = str2num(split{j});
        end

    end
    
    % Close file for reading
    fclose(fid);
end


