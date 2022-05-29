function [mesh] = gen2DMesh(filename)
    % Generates 2D mesh information reading
    % from GMSH mesh format

    % (C) 2016 - 2017 - Juan Marcelo Gimenez
    % (C) 2017 - Santiago Marquez Damian
    %            - Formatting
    %            - Native GMSH data support
    %            - Element to face map added (mesh.ele)

    % Raw data
    %mesh.xnod = load('xnod2');
    %mesh.icone = load('icone2');

    % GMSH reading
    [mesh.xnod, mesh.icone] = readGMSH(filename);
    % Override third dimension in points
    mesh.xnod = mesh.xnod(:, 1:2);

    [mesh.nelem, mesh.nen] = size(mesh.icone);
    [mesh.nnod, mesh.ndim] = size(mesh.xnod);

    % Centroids
    mesh.C = zeros(mesh.nelem,2);
    for i = 1:mesh.nen
        mesh.C(:,1) = mesh.C(:,1) + mesh.xnod(mesh.icone(:,i),1);
        mesh.C(:,2) = mesh.C(:,2) + mesh.xnod(mesh.icone(:,i),2);
    end
    mesh.C = mesh.C./mesh.nen;

    % figure(1);clf;pltmsh(mesh.xnod,mesh.icone);axis equal;
    % hold on; plot(mesh.C(:,1),mesh.C(:,2),'*g');

    mesh.faces = [];
    mesh.owner = [];
    mesh.neighbour = [];
    mesh.Sf = [];
    mesh.Cf = [];
    delta = [];
    k = [];
    mesh.V = [];
    
    for i = 1:mesh.nelem
        v = 0;
        for j = 1:mesh.nen
            f = [mesh.icone(i,j), mesh.icone(i, mod(j,mesh.nen) + 1)];
            
            % Check if this face already exists
            mesh.iface = [];
            if mesh.faces
                mesh.iface = find(mesh.faces(:, 1) == f(2) & ...
                mesh.faces(:, 2) == f(1));
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
end
