function [A, b] = assemble_advection_hr(Mesh, patches, v, psif)

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

            if dot(v(iface, :),Mesh.Sf(iface, :)) > 0
              fx = norm(Mesh.C(n, :) - Mesh.Cf(iface, :),2)/d;
              A(o, o) = A(o, o) + (1 - 0.5*psif(iface))*dot(v(iface, :),Mesh.Sf(iface, :));
              A(o, n) = A(o, n) + 0.5*psif(iface)*dot(v(iface, :),Mesh.Sf(iface, :));
              A(n, o) = A(n, o) - (1 - 0.5*psif(iface))*dot(v(iface, :),Mesh.Sf(iface, :));
              A(n, n) = A(n, n) - 0.5*psif(iface)*dot(v(iface, :),Mesh.Sf(iface, :));
            else
              fx = norm(Mesh.C(o, :) - Mesh.Cf(iface, :),2)/d;
              A(n, n) = A(n, n) - (1 - 0.5*psif(iface))*dot(v(iface, :),Mesh.Sf(iface, :));
              A(n, o) = A(n, o) - 0.5*psif(iface)*dot(v(iface, :),Mesh.Sf(iface, :));
              A(o, n) = A(o, n) + (1 - 0.5*psif(iface))*dot(v(iface, :),Mesh.Sf(iface, :));
              A(o, o) = A(o, o) + 0.5*psif(iface)*dot(v(iface, :),Mesh.Sf(iface, :));
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
              if dot(v(iface, :),Mesh.Sf(iface, :)) > 0
                b(o) = b(o) - phib.value*dot(v(iface, :),Mesh.Sf(iface, :));
              else
                b(o) = b(o) - phib.value*dot(v(iface, :),Mesh.Sf(iface, :));
              end
            elseif strcmp(phib.type, 'Neumann')
              d = norm(Mesh.C(o, :) - Mesh.Cf(iface, :),2);
              if dot(v(iface, :),Mesh.Sf(iface, :)) > 0

                A(o, o) = A(o, o) + dot(v(iface, :),Mesh.Sf(iface, :));
                b(o) = b(o) - phib.value*d*dot(v(iface, :),Mesh.Sf(iface, :));
              else
                A(o, o) = A(o, o) + dot(v(iface, :),Mesh.Sf(iface, :));
                b(o) = b(o) + phib.value*d*dot(v(iface, :),Mesh.Sf(iface, :));
              end
            end
         end
     end
end
