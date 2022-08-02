function [sum_flux_cell] = sum_flux_cells(flux, Mesh)
    sum_flux_cell = zeros(Mesh.ncells, 1);
    for iface = 1:Mesh.nfaces
         o = Mesh.owner(iface);
         n = Mesh.neighbour(iface);
         if n>0
            sum_flux_cell(o, :) = sum_flux_cell(o, :) + flux(iface);
            sum_flux_cell(n, :) = sum_flux_cell(n, :) - flux(iface);
         else
            sum_flux_cell(o, :) = sum_flux_cell(o, :) + flux(iface);
         end
     endfor
end
