function b = assemble_diffusion(Mesh, Q)

    b = zeros(Mesh.ncells, 1);

    if (length(Q) == 1)
        Q = Q*ones(Mesh.nelem, 1);
    end
    for i = 1:Mesh.nelem
      b(i) = Mesh.V(i)*Q(i);
    end
end
