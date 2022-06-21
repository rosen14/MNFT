function Mesh = interp_phif(Mesh, patches, phi_k, v, tvd_sch)
    Mesh.phif = zeros(Mesh.nfaces, 1);
    Mesh.rf = zeros(Mesh.nfaces, 1);
    Mesh.psif = zeros(Mesh.nfaces, 1);
    Mesh.grad = zeros(Mesh.ncells, 2);
    for iface = 1:Mesh.nfaces
         o = Mesh.owner(iface);
         n = Mesh.neighbour(iface);
         S = norm(Mesh.Sf(iface, :), 2);
         if n>0
            % tratamiento de caras internas
            dv = Mesh.C(n, :) - Mesh.C(o, :);
            d = norm(dv, 2);
            fx = norm(Mesh.C(o, :) - Mesh.Cf(iface, :),2)/d;
            Mesh.phif(iface) = phi_k(n)*fx + (1-fx)*phi_k(o);

            Mesh.grad(o, :) = Mesh.grad(o, :) + Mesh.phif(iface)*Mesh.Sf(iface, :)/Mesh.V(o);
            Mesh.grad(n, :) = Mesh.grad(n, :) - Mesh.phif(iface)*Mesh.Sf(iface, :)/Mesh.V(n);
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
                Mesh.phif(iface) = phib.value;
                Mesh.grad(o, :) = Mesh.grad(o, :) + Mesh.phif(iface)*Mesh.Sf(iface, :)/Mesh.V(o);

            end
            if strcmp(phib.type, 'Neumann')
                d_of = norm(Mesh.C(o, :) - Mesh.Cf(iface, :),2);
                Mesh.phif(iface) =  phi_k(o) + phib.value*d_of;
                Mesh.grad(o, :) = Mesh.grad(o, :) + Mesh.phif(iface)*Mesh.Sf(iface, :)/Mesh.V(o);
            end
            if strcmp(phib.type, 'Robin')
                %completar aqui
            end
         end
     endfor

    for iface = 1:Mesh.nfaces
         o = Mesh.owner(iface);
         n = Mesh.neighbour(iface);
         S = norm(Mesh.Sf(iface, :), 2);
         if n>0
           if dot(v(iface, :),Mesh.Sf(iface, :)) > 0
             Mesh.rf(iface) = 2*dot(Mesh.grad(o, :), dv)/(phi_k(n) - phi_k(o) + 1e-12) - 1;
           else
             Mesh.rf(iface) = 2*dot(Mesh.grad(n, :), -dv)/(phi_k(o) - phi_k(n) + 1e-12) - 1;
           end
           if strcmp(tvd_sch, 'superbee')
             Mesh.psif(iface) = max(0, max(min(1, 2*Mesh.rf(iface)), min(2, Mesh.rf(iface))));
           elseif strcmp(tvd_sch, 'minmod')
             Mesh.psif(iface) = max(0, min(1, Mesh.rf(iface)));
           elseif strcmp(tvd_sch, 'upwind')
             Mesh.psif(iface) = 0;
           elseif strcmp(tvd_sch, 'central_diff')
             Mesh.psif(iface) = 1;
           endif
         endif
     end
end
