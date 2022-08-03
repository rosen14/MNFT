function [A, b] = assemble_advection_hr_flux(Mesh, patches, flux, psif)

    A = sparse(zeros(Mesh.ncells, Mesh.ncells));
    b = zeros(Mesh.ncells, 1);

    for iface = 1:Mesh.nfaces
         o = Mesh.owner(iface);
         n = Mesh.neighbour(iface);
         S = norm(Mesh.Sf(iface, :), 2);
         if n>0
            % tratamiento de caras internas
            dv = Mesh.C(n, :) - Mesh.C(o, :);
            d = norm(dv, 2);

            if flux(iface) > 0
              fx = norm(Mesh.C(n, :) - Mesh.Cf(iface, :),2)/d;
              A(o, o) = A(o, o) + (1 - 0.5*psif(iface))*flux(iface);
              A(o, n) = A(o, n) + 0.5*psif(iface)*flux(iface);
              A(n, o) = A(n, o) - (1 - 0.5*psif(iface))*flux(iface);
              A(n, n) = A(n, n) - 0.5*psif(iface)*flux(iface);
            else
              fx = norm(Mesh.C(o, :) - Mesh.Cf(iface, :),2)/d;
              A(n, n) = A(n, n) - (1 - 0.5*psif(iface))*flux(iface);
              A(n, o) = A(n, o) - 0.5*psif(iface)*flux(iface);
              A(o, n) = A(o, n) + (1 - 0.5*psif(iface))*flux(iface);
              A(o, o) = A(o, o) + 0.5*psif(iface)*flux(iface);
            endif


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
              if flux(iface) > 0
                b(o) = b(o) - phib.value*flux(iface);
              else
                b(o) = b(o) - phib.value*flux(iface);
              end
            elseif strcmp(phib.type, 'Neumann')
              d = norm(Mesh.C(o, :) - Mesh.Cf(iface, :),2);
              if flux(iface) > 0

                A(o, o) = A(o, o) + flux(iface);
                b(o) = b(o) - phib.value*d*flux(iface);
              else
                A(o, o) = A(o, o) + flux(iface);
                b(o) = b(o) + phib.value*d*flux(iface);
              end
            end
         end
     end
end
