function [A, b] = assemble_advection(Mesh, patches, v, aprox)

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

            if aprox == 'UW'
              A(o, o) = A(o, o) + dot(v(iface, :),Mesh.Sf(iface, :));
              A(n, o) = A(n, o) - dot(v(iface, :),Mesh.Sf(iface, :));
            elseif aprox == 'CD'
              fx = norm(Mesh.C(n, :) - Mesh.Cf(iface, :),2)/d;
              A(o, o) = A(o, o) + fx*dot(v(iface, :),Mesh.Sf(iface, :));
              A(o, n) = A(o, n) + (1-fx)*dot(v(iface, :),Mesh.Sf(iface, :));
              A(n, o) = A(n, o) - fx*dot(v(iface, :),Mesh.Sf(iface, :));
              A(n, n) = A(n, n) - (1-fx)*dot(v(iface, :),Mesh.Sf(iface, :));
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
            if aprox == 'UW'
              if dot(v(iface, :),Mesh.Sf(iface, :)) > 0
                A(o, o) = A(o, o) + dot(v(iface, :),Mesh.Sf(iface, :));
              else
                b(o) = b(o) + phib.value*dot(v(iface, :),Mesh.Sf(iface, :));
              end
            end
            if aprox == 'CD'
              if dot(v(iface, :),Mesh.Sf(iface, :)) > 0
                b(o) = b(o) - phib.value*dot(v(iface, :),Mesh.Sf(iface, :));
              else
                b(o) = b(o) + phib.value*dot(v(iface, :),Mesh.Sf(iface, :));
              end
             end
         end
end
