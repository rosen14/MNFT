function [psif] = get_psif(Mesh, patches, phi_k, v, phif, grad, tvd_sch)

    rf = zeros(Mesh.nfaces, 1);
    psif = zeros(Mesh.nfaces, 1);

    for iface = 1:Mesh.nfaces
         o = Mesh.owner(iface);
         n = Mesh.neighbour(iface);
         S = norm(Mesh.Sf(iface, :), 2);
         if n>0
           dv = Mesh.C(n, :) - Mesh.C(o, :);
           if dot(v(iface, :),Mesh.Sf(iface, :)) > 0
             rf(iface) = 2*dot(grad(o, :), dv)/(phi_k(n) - phi_k(o) + 1e-12) - 1;
           else
             rf(iface) = 2*dot(grad(n, :), -dv)/(phi_k(o) - phi_k(n) + 1e-12) - 1;
           end
           if strcmp(tvd_sch, 'superbee')
             psif(iface) = max(0, max(min(1, 2*rf(iface)), min(2, rf(iface))));
           elseif strcmp(tvd_sch, 'minmod')
             psif(iface) = max(0, min(1, rf(iface)));
           elseif strcmp(tvd_sch, 'upwind')
             psif(iface) = 0;
           elseif strcmp(tvd_sch, 'central_diff')
             psif(iface) = 1;
           endif
         endif
     end
end
