function [A, b] = assemble_diffusion(Mesh, patches, nu)

    A = sparse(zeros(Mesh.ncells, Mesh.ncells));
    b = zeros(Mesh.ncells, 1);
    
    if (length(nu) == 1)
        nu = nu*ones(Mesh.nfaces, 1);
    end
    for iface = 1:Mesh.nfaces
         o = Mesh.owner(iface);
         n = Mesh.neighbour(iface);
         S = norm(Mesh.Sf(iface, :), 2);
         if n>0
            % tratamiento de caras internas
            dv = Mesh.C(n, :) - Mesh.C(o, :);
            d = norm(dv, 2);

            A(o, o) = A(o, o) + nu(iface)/d*S;
            A(o, n) = A(o, n) - nu(iface)/d*S;
            A(n, o) = A(n, o) - nu(iface)/d*S;
            A(n, n) = A(n, n) + nu(iface)/d*S;
         else
            % tratamiento de caras de frontera        
            phib = -1;
            for ipatch = 1:length(patches)
                if find(patches{ipatch}.faces == iface) 
                    phib = patches{ipatch};
                    break;
                end
            end
            if strcmp(phib.type, 'Dirichlet')
                dv = Mesh.Cf(iface, :) - Mesh.C(o, :);
                dn = dv*Mesh.Sf(iface, :)'/S*(Mesh.Sf(iface, :)/S);
                d = norm(dn, 2); 
                A(o, o) = A(o, o) + nu(iface)/d*S;
                b(o) = b(o) + phib.value*nu(iface)/d*S;
            end
            if strcmp(phib.type, 'Neumann')
                b(o) = b(o) + phib.value;
            end
            if strcmp(phib.type, 'Robin')
                %completar aqui
            end
         end
     end
end
