function [divphi] = fvc_div(phif, Mesh)
    divphi = zeros(Mesh.ncells, 1);
    for iface = 1:Mesh.nfaces
         o = Mesh.owner(iface);
         n = Mesh.neighbour(iface);
         if n>0
            divphi(o, :) = divphi(o, :) + dot(phif(iface,:),Mesh.Sf(iface, :));
            divphi(n, :) = divphi(n, :) - dot(phif(iface,:),Mesh.Sf(iface, :));
         else
            divphi(o, :) = divphi(o, :) + dot(phif(iface,:),Mesh.Sf(iface, :));
         end
     endfor
end
