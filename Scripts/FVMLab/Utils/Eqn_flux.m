function [flux] = Eqn_flux(Mesh, patches, nu, A, phi)
   flux = zeros(Mesh.nfaces, 1);
   for iface = 1:Mesh.nfaces
         o = Mesh.owner(iface);
         n = Mesh.neighbour(iface);
         if n>0
           flux(iface) = A(o, n)*(phi(n) - phi(o));
         endif
     end
endfunction
