function [phif, grad] = interp_phif_v2(Mesh, patches, phi_k)
    phif = zeros(Mesh.nfaces, 1);
    grad = zeros(Mesh.ncells, 2);
    for iface = 1:Mesh.nfaces
         o = Mesh.owner(iface);
         n = Mesh.neighbour(iface);
         S = norm(Mesh.Sf(iface, :), 2);
         if n>0
            % Calculo phif siguiendo el gf de pag. 160 Moukalled
            % tratamiento de caras internas
            ef = Mesh.Sf(iface,:)/norm(Mesh.Sf(iface,:),2); % surface unit vector
            d_Cf = Mesh.Cf(iface, :) - Mesh.C(o, :);
            d_fF = Mesh.C(n, :) - Mesh.Cf(iface, :);
            dv = Mesh.C(n, :) - Mesh.C(o, :);
            d = norm(dv, 2);
            gf = dot(d_Cf, ef)/(dot(d_Cf, ef) + dot(d_fF, ef));

            phif(iface) = phi_k(n)*gf + (1-gf)*phi_k(o);

            grad(o, :) = grad(o, :) + phif(iface)*Mesh.Sf(iface, :)/Mesh.V(o);
            grad(n, :) = grad(n, :) - phif(iface)*Mesh.Sf(iface, :)/Mesh.V(n);
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
                phif(iface) = phib.value;
                grad(o, :) = grad(o, :) + phif(iface)*Mesh.Sf(iface, :)/Mesh.V(o);

            end
            if strcmp(phib.type, 'Neumann')
                d_of = norm(Mesh.C(o, :) - Mesh.Cf(iface, :),2);
                phif(iface) =  phi_k(o) + phib.value*d_of;
                grad(o, :) = grad(o, :) + phif(iface)*Mesh.Sf(iface, :)/Mesh.V(o);
            end
            if strcmp(phib.type, 'Robin')
                %completar aqui
            end
         end
     endfor
end
