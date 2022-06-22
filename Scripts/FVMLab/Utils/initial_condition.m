function phi_init = initial_condition(Mesh)

   phi_init = zeros(Mesh.ncells,1);

   for cell = 1:Mesh.ncells
      if ((1/6 < Mesh.C(cell, 1)) && (Mesh.C(cell, 1) < 0.5))
        phi_init(cell) = sin(3*pi*(Mesh.C(cell, 1)-1/6))**2;
      endif
   endfor
