function [flux] = Eqn_flux(Mesh, patches, nu, A, phi)
   flux = zeros(Mesh.nfaces, 1);
   for iface = 1:Mesh.nfaces
         o = Mesh.owner(iface);
         n = Mesh.neighbour(iface);
         S = norm(Mesh.Sf(iface, :), 2);
         if n>0
           flux(iface) = A(o, n)*(phi(n) - phi(o));
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
                flux(iface) = nu(iface)*S*(phib.value - phi(o))/d;
            end
            if strcmp(phib.type, 'Neumann')
                flux(iface) = nu(iface)*S*phib.value;
            end
            if strcmp(phib.type, 'Robin')
                %completar aqui
            end
         end
     end
endfunction
